process PREPARE_FLANKS {
    tag "${meta.id}"
    label 'python'
    publishDir "${params.outdir}/genomes/${meta.id}", mode: 'copy', pattern: 'flanks.*'
    publishDir "${params.outdir}/genomes/${meta.id}", mode: 'copy', pattern: 'targets.tsv'
    input:
    tuple val(meta), path(fasta, stageAs: 'genome.fa'), path(fai, stageAs: 'genome.fa.fai'), path(annotation, stageAs: 'annotation.input'), val(format), path(chrom_map, stageAs: 'chrom_map.json')
    path code_files, stageAs: 'code/*'
    val options
    output:
    tuple val(meta), path('flanks.fa'), path('flanks.tsv'), emit: flanks
    tuple val(meta), path('targets.tsv'), emit: targets
    path 'flanks.python_packages.txt', emit: versions
    script:
    def q = { x -> "'" + x.toString().replace("'", "'\"'\"'") + "'" }
    def inputArg = format == 'rm_out' ? '--rm-out annotation.input' : "--bed annotation.input --bed-format ${q(format)}"
    def mapArg = chrom_map ? '--chrom-map chrom_map.json' : ''
    def fields = ['flank_len','nrt_len','overhang_len','target_near_keep','repeat_end_keep','max_groups','target_groups','target_repeats']
    def args = fields.collect { "--${it.replace('_', '-')} ${q(options[it])}" }.join(' ')
    if (options.soft_mask) args += ' --soft-mask'
    """
    python code/prepare_flanks.py ${inputArg} --genome-fasta genome.fa ${mapArg} ${args} --out-table flanks.tsv --out-fasta flanks.fa --out-targets targets.tsv
    python -c "import importlib.metadata as m; print('\\n'.join(sorted(str(d.metadata['Name']) + '==' + d.version for d in m.distributions())))" > flanks.python_packages.txt
    """
}
