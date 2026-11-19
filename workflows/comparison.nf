include { REPEATMASKER } from '../modules/repeatmasker'
include { FASTA_INDEX } from '../modules/fasta_index'
include { PREPARE_REPEAT_ANNOTATION } from '../modules/prepare_repeat_annotation'
include { PREPARE_FLANKS } from '../modules/prepare_flanks'
include { MINIMAP_INDEX } from '../modules/minimap_index'
include { MINIMAP_ALIGN } from '../modules/minimap_align'
include { CLASSIFY_ALIGNMENTS } from '../modules/classify_alignments'
include { VALIDATE_HAL } from '../modules/validate_hal'

// Manifests use literal TSV fields without CSV quoting or embedded tabs/newlines.
def readManifest(input, required) {
    if (!input) error 'Both --genomes and --pairs manifests are required; see README.md'
    def path = file(input, checkIfExists: true)
    def lines = path.readLines().findAll { it.trim() && !it.startsWith('#') }
    if (!lines) error "Empty manifest: ${path}"
    def columns = lines[0].split('\t', -1).toList()
    if (columns.unique(false).size() != columns.size()) error "Duplicate manifest headers: ${path}"
    required.each { if (!columns.contains(it)) error "Missing ${it} column in ${path}" }
    def rows = lines.drop(1).collect { line ->
        def fields = line.split('\t', -1).toList()
        if (fields.size() != columns.size()) error "Wrong number of TSV fields in ${path}"
        def row = [:]
        columns.eachWithIndex { name, i -> row[name] = fields[i].trim() }
        row
    }
    if (!rows) error "No data rows in ${path}"
    [path, rows]
}

def resolveInput(base, value) {
    if (!value) return []
    def p = java.nio.file.Paths.get(value)
    file((p.isAbsolute() ? p : base.resolve(p)).normalize().toString(), checkIfExists: true)
}

def requireInput(value, message) {
    if (!value) error message
}

def checkArgs(values, blocked, label) {
    if (!(values instanceof List) || values.any { !(it instanceof String) || it.contains('\n') }) {
        error "${label} must be a YAML/JSON list of argument strings"
    }
    values.each { token ->
        if (blocked.any { token == it || token.startsWith(it + '=') || (it.size() == 2 && token.startsWith(it)) }) {
            error "Reserved option ${token} in ${label}; use the named configuration parameter"
        }
    }
}

workflow COMPARISON {
    main:
    def starts = [genome:0, annotation:1, flanks:2, paf:3, statuses:4]
    def ends = [annotation:1, flanks:2, paf:3, statuses:4, validation:5]
    if (!starts.containsKey(params.from) || !ends.containsKey(params.until)) error 'Invalid --from or --until'
    def first = starts[params.from]
    def last = ends[params.until]
    if (last < first || (first >= 2 && last == first)) error 'Requested stage range has no executable work'
    if (!(params.rm_sensitivity in ['default','s','q','qq'])) error 'rm_sensitivity must be default, s, q or qq'
    if (!(params.minimap_secondary in ['yes','no'])) error 'minimap_secondary must be yes or no'
    if (!(params.classify_mode in ['stream','collect'])) error 'classify_mode must be stream or collect'
    if (params.rm_engine == 'hmmer' && params.rm_cpus < 2) error 'HMMER requires at least two allocated CPUs per batch'
    checkArgs(params.minimap_index_args, ['-a','-o','-d','-x','-t','--split-prefix'], 'minimap_index_args')
    checkArgs(params.minimap_extra_args, ['-a','-o','-d','-x','-t','-N','-p','-c','-k','-w','-H','-I','--cs','--secondary','--split-prefix'], 'minimap_extra_args')
    checkArgs(params.rm_extra_args, ['-engine','-species','-lib','-pa','-dir','-s','-q','-qq','-align'], 'rm_extra_args')

    def genomeInput = readManifest(params.genomes, ['genome_id'])
    def pairInput = readManifest(params.pairs, ['source','target'])
    def genomePath = genomeInput[0]
    def pairPath = pairInput[0]
    def genomes = genomeInput[1].collect { row ->
        if (!(row.genome_id ==~ /[A-Za-z][A-Za-z0-9_]*/)) error "Invalid genome_id: ${row.genome_id}"
        def g = [id:row.genome_id, species:row.species ?: '', format:row.annotation_format ?: '']
        ['fasta','annotation','flank_fasta','flank_table','chrom_map'].each {
            g[it] = resolveInput(genomePath.parent, row[it])
        }
        if (g.annotation && !(g.format in ['rm_out','bed6','custom5'])) error "Set annotation_format for ${g.id}"
        g
    }
    if (genomes.collect { it.id }.unique(false).size() != genomes.size()) error 'Duplicate genome_id'
    def genomeMap = genomes.collectEntries { [(it.id):it] }
    def pairs = pairInput[1].collect { row ->
        if (!genomeMap.containsKey(row.source) || !genomeMap.containsKey(row.target)) error 'Unknown source or target genome'
        if (row.source == row.target) error 'Self comparisons are not supported in the pair manifest'
        def p = [source:row.source, target:row.target, preset:row.preset ?: params.minimap_preset,
                 hal_source:row.hal_source ?: '', hal_target:row.hal_target ?: '']
        if (!(p.preset ==~ /[A-Za-z0-9_-]+/)) error "Invalid preset: ${p.preset}"
        ['paf','statuses','hal','hal_source_chrom_map','hal_target_chrom_map'].each {
            p[it] = resolveInput(pairPath.parent, row[it])
        }
        p
    }
    if (pairs.collect { [it.source,it.target] }.unique(false).size() != pairs.size()) error 'Duplicate directed genome pair'
    if (pairs.collect { it.source + '__to__' + it.target }.unique(false).size() != pairs.size()) error 'Pair IDs collide in output directory names; choose unambiguous genome IDs'
    def sourceIds = pairs.collect { it.source }.unique()
    def sourceGenomes = sourceIds.collect { genomeMap[it] }
    def annotationIds = params.annotate_targets ? (sourceIds + pairs.collect { it.target }).unique() : sourceIds
    def annotationGenomes = annotationIds.collect { genomeMap[it] }
    if (params.annotate_targets && first > 1) error 'annotate_targets applies only to genome/annotation starts'
    if (first <= 1) {
        annotationGenomes.each { g ->
            if (!g.annotation) {
                if (first == 1) error "Missing external annotation for ${g.id}"
                requireInput(g.fasta, "FASTA required for RepeatMasker: ${g.id}")
                if (!params.rm_library) requireInput(g.species, "species required for RepeatMasker: ${g.id}")
            }
        }
        if (last >= 2) sourceGenomes.each { requireInput(it.fasta, "FASTA required for flanks: ${it.id}") }
    }
    if (first == 2) sourceGenomes.each { requireInput(it.flank_fasta, "flank_fasta required for ${it.id}") }
    if (first >= 2) sourceGenomes.each { requireInput(it.flank_table, "flank_table required for ${it.id}") }
    if (first <= 2 && last >= 3) pairs.each { requireInput(genomeMap[it.target].fasta, "Target FASTA required: ${it.target}") }
    if (first == 3) pairs.each { requireInput(it.paf, "External PAF required: ${it.source}->${it.target}") }
    if (first == 4) pairs.each { requireInput(it.statuses, "External statuses required: ${it.source}->${it.target}") }
    if (last == 5) pairs.each {
        requireInput(it.hal, 'HAL file required for every pair when --until validation')
        requireInput(it.hal_source, 'HAL source genome name required')
        requireInput(it.hal_target, 'HAL target genome name required')
    }
    genomes.findAll { it.fasta }.each {
        if (it.fasta.name.endsWith('.gz')) error "Use uncompressed FASTA: ${it.id}"
    }
    def library = params.rm_library ? file(params.rm_library, checkIfExists:true) : []
    def code = file("${projectDir}/bin/*.py", checkIfExists:true)

    // Keep unrelated parameters out of each process's task cache key.
    def flank_options = ['flank_len','nrt_len','overhang_len','target_near_keep','repeat_end_keep','max_groups','target_groups','target_repeats','soft_mask']
        .collectEntries { [(it):params[it]] }
    def classify_options = params.keySet().findAll { it.startsWith('classify_') }.sort()
        .collectEntries { [(it):params[it]] }
    def hal_options = ['hal_anchor_len','hal_anchor_offsets','hal_min_pred_recall','hal_alt_iou_threshold','hal_len_ratio_max','hal_batch_size','hal_liftover','hal_version']
        .collectEntries { [(it):params[it]] }

    // Record inputs and effective settings without reading genomic datasets.
    def stamp = new Date().format('yyyyMMdd-HHmmss-SSS')
    def provenance = file("${params.outdir}/provenance/${stamp}")
    java.nio.file.Files.createDirectories(provenance)
    java.nio.file.Files.copy(genomePath, provenance.resolve('genomes.tsv'))
    java.nio.file.Files.copy(pairPath, provenance.resolve('pairs.tsv'))
    provenance.resolve('parameters.json').text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(params))
    provenance.resolve('run.txt').text = "Nextflow: ${workflow.nextflow.version}\nCommand: ${workflow.commandLine}\nGenome manifest: ${genomePath}\nPair manifest: ${pairPath}\n"
    log.info "ComRepka: ${pairs.size()} pairs; ${params.from} -> ${params.until}"

    flank_ch = Channel.empty()
    if (first <= 1) {
        rm_ch = Channel.fromList(annotationGenomes.findAll { !it.annotation }.collect { tuple([id:it.id,species:it.species], it.fasta, library) })
        REPEATMASKER(rm_ch)
        annotation_ch = Channel.fromList(annotationGenomes.findAll { it.annotation }.collect { tuple([id:it.id],it.annotation,it.format,it.chrom_map) })
            .mix(REPEATMASKER.out.annotation.map { meta, annotation -> tuple(meta,annotation,'rm_out',genomeMap[meta.id].chrom_map) })
        PREPARE_REPEAT_ANNOTATION(annotation_ch, code)
        if (last >= 2) {
            FASTA_INDEX(Channel.fromList(sourceGenomes.collect { tuple([id:it.id],it.fasta) }))
            flank_inputs = FASTA_INDEX.out.indexed.map { meta, fai -> tuple(meta.id,meta,genomeMap[meta.id].fasta,fai) }
                .join(annotation_ch.filter { meta, annotation, format, cmap -> sourceIds.contains(meta.id) }
                    .map { meta, annotation, format, cmap -> tuple(meta.id,annotation,format,cmap) })
                .map { id,meta,fasta,fai,annotation,format,cmap -> tuple(meta,fasta,fai,annotation,format,cmap) }
            PREPARE_FLANKS(flank_inputs, code, flank_options)
            flank_ch = PREPARE_FLANKS.out.flanks.map { meta,fasta,table -> tuple(meta.id,fasta,table) }
        }
    } else {
        flank_ch = Channel.fromList(sourceGenomes.collect { tuple(it.id,it.flank_fasta,it.flank_table) })
    }

    paf_ch = Channel.empty()
    if (last >= 3 && first <= 2) {
        index_requests = pairs.collect { [target:it.target,preset:it.preset] }.unique().collect {
            tuple(it, genomeMap[it.target].fasta, params.minimap_index_args)
        }
        MINIMAP_INDEX(Channel.fromList(index_requests))
        alignment_inputs = Channel.fromList(pairs).map { tuple(it.source,it) }
            .combine(flank_ch, by:0)
            .map { source,meta,fasta,table -> tuple("${meta.target}:${meta.preset}",meta,fasta) }
            .combine(MINIMAP_INDEX.out.index.map { meta,index -> tuple("${meta.target}:${meta.preset}",index) }, by:0)
            .map { key,meta,fasta,index -> tuple(meta,index,fasta,params.minimap_index_args,params.minimap_extra_args) }
        MINIMAP_ALIGN(alignment_inputs)
        paf_ch = MINIMAP_ALIGN.out.paf
    } else if (first == 3) {
        paf_ch = Channel.fromList(pairs.collect { tuple(it,it.paf) })
    }
    status_ch = Channel.empty()
    if (last >= 4 && first <= 3) {
        classify_inputs = paf_ch.map { meta,paf -> tuple(meta.source,meta,paf) }
            .combine(flank_ch.map { source,fasta,table -> tuple(source,table) }, by:0)
            .map { source,meta,paf,table -> tuple(meta,paf,table) }
        CLASSIFY_ALIGNMENTS(classify_inputs, code, classify_options)
        status_ch = CLASSIFY_ALIGNMENTS.out.statuses
    } else if (first == 4) {
        status_ch = Channel.fromList(pairs.collect { tuple(it,it.statuses) })
    }
    if (last == 5) {
        validation_inputs = status_ch.map { meta,statuses -> tuple(meta.source,meta,statuses) }
            .combine(flank_ch.map { source,fasta,table -> tuple(source,table) }, by:0)
            .map { source,meta,statuses,table -> tuple(meta,statuses,table,meta.hal,meta.hal_source_chrom_map,meta.hal_target_chrom_map) }
        VALIDATE_HAL(validation_inputs, code, hal_options)
    }
}
