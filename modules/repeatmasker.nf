process REPEATMASKER {
    tag "${meta.id}"
    label 'repeatmasker'
    publishDir "${params.outdir}/genomes/${meta.id}", mode: 'copy', pattern: 'repeatmasker'
    input:
    tuple val(meta), path(fasta, stageAs: 'genome.fa'), path(library, stageAs: 'library.fa')
    output:
    tuple val(meta), path('repeatmasker/genome.fa.out'), emit: annotation
    path 'repeatmasker', emit: files
    script:
    def q = { x -> "'" + x.toString().replace("'", "'\"'\"'") + "'" }
    def pa = params.rm_engine == 'hmmer' ? Math.max(1, task.cpus.intdiv(2)) : task.cpus
    def sensitivity = params.rm_sensitivity == 'default' ? '' : "-${params.rm_sensitivity}"
    def reference = library ? "-lib ${q(library)}" : "-species ${q(meta.species)}"
    def extra = params.rm_extra_args.collect { q(it) }.join(' ')
    """
    mkdir repeatmasker
    RepeatMasker -engine ${q(params.rm_engine)} ${sensitivity} -align ${reference} -pa ${pa} -dir repeatmasker ${extra} genome.fa 2>&1 | tee repeatmasker/run.log
    RepeatMasker -v > repeatmasker/software_versions.txt 2>&1
    printf '%s\\n' ${q(params.rm_reference_id)} >> repeatmasker/software_versions.txt
    """
}
