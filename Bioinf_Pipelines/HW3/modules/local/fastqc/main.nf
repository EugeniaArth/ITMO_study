process FASTQC {
    tag "$sample_id ($stage)"
    label 'process_low'

    publishDir { "${params.outdir}/fastqc_${stage}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(stage), path(read1), path(read2)

    output:
    tuple val(sample_id), val(stage), path("*_fastqc.zip"), path("*_fastqc.html")

    script:
    """
    fastqc -t ${task.cpus} ${read1} ${read2}
    """
}
