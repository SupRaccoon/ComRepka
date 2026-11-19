process VALIDATE_HAL {
    tag "${meta.source}->${meta.target}"
    label 'hal'
    publishDir "${params.outdir}/validation/${meta.source}__to__${meta.target}", mode: 'copy'
    input:
    tuple val(meta), path(statuses, stageAs: 'statuses.tsv'), path(flank_table, stageAs: 'flanks.tsv'), path(hal, stageAs: 'alignment.hal'), path(source_map, stageAs: 'source_map.json'), path(target_map, stageAs: 'target_map.json')
    path code_files, stageAs: 'code/*'
    val options
    output:
    tuple val(meta), path('hal_statuses.tsv'), emit: statuses
    tuple path('hal_ambiguous.tsv'), path('hal_disagree.tsv'), path('hal_other_contigs.tsv'), path('hal_unsupported.tsv'), emit: diagnostics
    path 'hal_version.txt', emit: version
    script:
    def q = { x -> "'" + x.toString().replace("'", "'\"'\"'") + "'" }
    def maps = (source_map ? '--src-chrom-map source_map.json ' : '') + (target_map ? '--dst-chrom-map target_map.json' : '')
    def fields = ['anchor_len','anchor_offsets','min_pred_recall','alt_iou_threshold','len_ratio_max','batch_size']
    def args = fields.collect { "--${it.replace('_', '-')} ${q(options['hal_' + it])}" }.join(' ')
    """
    python code/validate_hal.py --flank-table flanks.tsv --pairs-tsv statuses.tsv --hal alignment.hal --hal-liftover ${q(options.hal_liftover)} --src-genome ${q(meta.hal_source)} --dst-genome ${q(meta.hal_target)} ${maps} ${args} --out-main hal_statuses.tsv --out-amb hal_ambiguous.tsv --out-dis hal_disagree.tsv --out-other hal_other_contigs.tsv --out-unsup hal_unsupported.tsv
    printf '%s\\n' ${q(options.hal_version)} > hal_version.txt
    """
}
