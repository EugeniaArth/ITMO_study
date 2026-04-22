nextflow.enable.dsl=2

params.reads   = "/Users/eugenianikonorova/Documents/ITMO/ITMO_study/data/*_{1,2}.fastq"
params.outdir  = "HW2/results"
params.threads = 4
params.memory  = "8 GB"

process FASTQC_RAW {
    tag "$sample_id"
    publishDir "${params.outdir}/02_fastqc_raw", mode: 'copy'

    cpus params.threads
    memory params.memory

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("*_fastqc.zip"), path("*_fastqc.html")

    script:
    """
    fastqc -t ${task.cpus} ${read1} ${read2}
    """
}

process TRIM_READS {
    tag "$sample_id"
    publishDir "${params.outdir}/03_trimmed_reads", mode: 'copy'

    cpus params.threads
    memory params.memory

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_1.fastq.gz"), path("${sample_id}_trimmed_2.fastq.gz")
    path("${sample_id}.fastp.json")
    path("${sample_id}.fastp.html")

    script:
    """
    fastp \
      -i ${read1} \
      -I ${read2} \
      -o ${sample_id}_trimmed_1.fastq.gz \
      -O ${sample_id}_trimmed_2.fastq.gz \
      --thread ${task.cpus} \
      --json ${sample_id}.fastp.json \
      --html ${sample_id}.fastp.html
    """
}

process FASTQC_TRIMMED {
    tag "$sample_id"
    publishDir "${params.outdir}/04_fastqc_trimmed", mode: 'copy'

    cpus params.threads
    memory params.memory

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("*_fastqc.zip"), path("*_fastqc.html")

    script:
    """
    fastqc -t ${task.cpus} ${read1} ${read2}
    """
}

workflow QC_AND_TRIM {
    take:
    reads_ch

    main:
    raw_qc_ch = FASTQC_RAW(reads_ch)

    trimmed_result = TRIM_READS(reads_ch)
    trimmed_reads_ch = trimmed_result[0]

    trimmed_qc_ch = FASTQC_TRIMMED(trimmed_reads_ch)

    emit:
    raw_qc_ch
    trimmed_reads_ch
    trimmed_qc_ch
}

workflow {
    reads_ch = Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .map { sid, reads -> tuple(sid, reads[0], reads[1]) }

    QC_AND_TRIM(reads_ch)
}