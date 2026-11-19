process PREPARE_REPEAT_ANNOTATION {
    tag "${meta.id}"
    label 'python'
    publishDir "${params.outdir}/genomes/${meta.id}", mode: 'copy'
    input:
    tuple val(meta), path(annotation, stageAs: 'annotation.input'), val(format), path(chrom_map, stageAs: 'chrom_map.json')
    path code_files, stageAs: 'code/*'
    output:
    tuple val(meta), path('repeats.bed'), emit: bed
    path 'annotation.python_packages.txt', emit: versions
    script:
    def q = { x -> "'" + x.toString().replace("'", "'\"'\"'") + "'" }
    def mapArg = chrom_map ? '--chrom-map chrom_map.json' : ''
    """
    python code/prepare_repeat_annotation.py --input annotation.input --format ${q(format)} --output repeats.bed --target-groups ${q(params.target_groups)} --target-repeats ${q(params.target_repeats)} ${mapArg}
    python -c "import importlib.metadata as m; print('\\n'.join(sorted(str(d.metadata['Name']) + '==' + d.version for d in m.distributions())))" > annotation.python_packages.txt
    """
}
