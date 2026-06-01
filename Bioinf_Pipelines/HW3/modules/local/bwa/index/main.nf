process INDEX_REFERENCE {
    tag "${ref.baseName}"
    label 'process_low'

    publishDir "${params.outdir}/reference_index", mode: 'copy'

    input:
    path ref

    output:
    tuple path("reference.fa"), path("reference.fa.*")

    script:
    """
    cp ${ref} reference.fa
    bwa index reference.fa
    samtools faidx reference.fa
    """
}
