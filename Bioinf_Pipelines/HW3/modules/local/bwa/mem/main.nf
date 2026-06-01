process MAP_READS {
    tag "$sample_id"
    label 'process_high'

    publishDir "${params.outdir}/mapping", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)
    tuple path(ref), path(ref_index_files)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai")

    script:
    """
    bwa mem -t ${task.cpus} ${ref} ${read1} ${read2} \\
      | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam

    samtools index ${sample_id}.sorted.bam
    """
}
