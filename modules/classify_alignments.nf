process CLASSIFY_ALIGNMENTS {
    tag "${meta.source}->${meta.target}"
    label 'python'
    publishDir "${params.outdir}/pairs/${meta.source}__to__${meta.target}", mode: 'copy'
    input:
    tuple val(meta), path(paf, stageAs: 'alignments.paf'), path(flank_table, stageAs: 'flanks.tsv')
    path code_files, stageAs: 'code/*'
    val options
    output:
    tuple val(meta), path('statuses.tsv'), emit: statuses
    path 'python_packages.txt', emit: versions
    script:
    def q = { x -> "'" + x.toString().replace("'", "'\"'\"'") + "'" }
    def fields = options.keySet().sort()
    def args = fields.collect {
        def flag = '--' + it.substring(9).replace('_', '-')
        options[it] instanceof Boolean ? (options[it] ? flag : '') : "${flag} ${q(options[it])}"
    }.join(' ')
    """
    python code/classify_alignments.py --paf alignments.paf --flank-table flanks.tsv --out-table statuses.tsv ${args}
    python -c "import importlib.metadata as m; print('\\n'.join(sorted(str(d.metadata['Name']) + '==' + d.version for d in m.distributions())))" > python_packages.txt
    """
}
