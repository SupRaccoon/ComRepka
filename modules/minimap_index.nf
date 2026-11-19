process MINIMAP_INDEX {
    tag "${meta.target}:${meta.preset}"
    label 'minimap'
    publishDir "${params.outdir}/genomes/${meta.target}/indices/${meta.preset}", mode: 'copy', pattern: '*.txt'
    input:
    tuple val(meta), path(fasta, stageAs: 'target.fa'), val(index_args)
    output:
    tuple val(meta), path('target.mmi'), emit: index
    path 'minimap2.version.txt', emit: version
    script:
    def q = { x -> "'" + x.toString().replace("'", "'\"'\"'") + "'" }
    def extra = index_args.collect { q(it) }.join(' ')
    """
    minimap2 -x ${q(meta.preset)} -t ${task.cpus} ${extra} -d target.mmi target.fa
    minimap2 --version > minimap2.version.txt
    """
}
