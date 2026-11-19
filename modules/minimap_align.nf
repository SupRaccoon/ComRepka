process MINIMAP_ALIGN {
    tag "${meta.source}->${meta.target}"
    label 'minimap'
    publishDir "${params.outdir}/pairs/${meta.source}__to__${meta.target}", mode: 'copy'
    input:
    tuple val(meta), path(index, stageAs: 'target.mmi'), path(flanks, stageAs: 'flanks.fa'), val(index_args), val(mapping_args)
    output:
    tuple val(meta), path('alignments.paf'), emit: paf
    path 'minimap2.version.txt', emit: version
    script:
    def q = { x -> "'" + x.toString().replace("'", "'\"'\"'") + "'" }
    def extra = (index_args + mapping_args).collect { q(it) }.join(' ')
    """
    minimap2 -x ${q(meta.preset)} -t ${task.cpus} -c --cs=long --secondary=${q(params.minimap_secondary)} -N ${q(params.minimap_N)} -p ${q(params.minimap_p)} ${extra} target.mmi flanks.fa > alignments.paf
    minimap2 --version > minimap2.version.txt
    """
}
