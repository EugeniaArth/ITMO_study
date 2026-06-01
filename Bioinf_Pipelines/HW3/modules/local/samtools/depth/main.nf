process CALCULATE_COVERAGE {
    tag "$sample_id"
    label 'process_low'

    publishDir "${params.outdir}/coverage", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.coverage.txt")

    script:
    """
    samtools depth ${bam} > ${sample_id}.coverage.txt
    """
}
