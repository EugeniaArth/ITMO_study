process INDEX_REFERENCE {
    tag "${ref.baseName}"
    label 'process_low'

    publishDir "${params.outdir}/reference_index", mode: 'copy'

    input:
    path ref, stageAs: 'reference.fa'

    output:
    tuple path("reference.fa"), path("reference.fa.*")

    script:
    """
    bwa index reference.fa
    samtools faidx reference.fa
    """
}
