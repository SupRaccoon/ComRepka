process FASTA_INDEX {
    tag "${meta.id}"
    label 'python'
    publishDir "${params.outdir}/genomes/${meta.id}", mode: 'copy', pattern: '*.fai'
    publishDir "${params.outdir}/genomes/${meta.id}", mode: 'copy', pattern: 'samtools.version.txt'
    input:
    tuple val(meta), path(fasta, stageAs: 'genome.fa')
    output:
    tuple val(meta), path('genome.fa.fai'), emit: indexed
    path 'samtools.version.txt', emit: version
    script:
    """
    samtools faidx genome.fa
    samtools --version > samtools.version.txt
    """
}
