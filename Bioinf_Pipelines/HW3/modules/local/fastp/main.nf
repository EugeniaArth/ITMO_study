process TRIM_READS {
    tag "$sample_id"
    label 'process_medium'

    publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_1.fastq.gz"), path("${sample_id}_trimmed_2.fastq.gz"), emit: reads
    path("${sample_id}.fastp.json"), emit: json
    path("${sample_id}.fastp.html"), emit: html

    script:
    """
    fastp \\
      -i ${read1} \\
      -I ${read2} \\
      -o ${sample_id}_trimmed_1.fastq.gz \\
      -O ${sample_id}_trimmed_2.fastq.gz \\
      --thread ${task.cpus} \\
      --json ${sample_id}.fastp.json \\
      --html ${sample_id}.fastp.html
    """
}
