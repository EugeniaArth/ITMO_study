nextflow.enable.dsl=2

params {
    reads: String = "/Users/eugenianikonorova/Documents/ITMO/data/*_{1,2}.fastq"
    reference: Path = "/Users/eugenianikonorova/Documents/ITMO/data/NZ_CP012026.1.fasta"
    outdir: Path = "HW2/results"
    threads: Integer = 10
    memory: String = "10 GB"
}

process FASTQC {
    tag "$sample_id ($stage)"
    publishDir { "${params.outdir}/fastqc_${stage}" }, mode: 'copy'

    cpus params.threads
    memory params.memory
    conda "${projectDir}/envs/fastqc.yml"

    input:
    tuple val(sample_id), val(stage), path(read1), path(read2)

    output:
    tuple val(sample_id), val(stage), path("*_fastqc.zip"), path("*_fastqc.html")

    script:
    """
    fastqc -t ${task.cpus} ${read1} ${read2}
    """
}

process TRIM_READS {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

    cpus params.threads
    memory params.memory
    conda "${projectDir}/envs/fastp.yml"

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_1.fastq.gz"), path("${sample_id}_trimmed_2.fastq.gz"), emit: reads
    path("${sample_id}.fastp.json"), emit: json
    path("${sample_id}.fastp.html"), emit: html

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

process INDEX_REFERENCE {
    tag "${ref.baseName}"
    publishDir "${params.outdir}/reference_index", mode: 'copy'

    cpus 1
    memory "2 GB"
    conda "${projectDir}/envs/bwa_samtools.yml"

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

process MAP_READS {
    tag "$sample_id"
    publishDir "${params.outdir}/mapping", mode: 'copy'

    cpus params.threads
    memory params.memory
    conda "${projectDir}/envs/bwa_samtools.yml"

    input:
    tuple val(sample_id), path(read1), path(read2)
    tuple path(ref), path(ref_index_files)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai")

    script:
    """
    bwa mem -t ${task.cpus} ${ref} ${read1} ${read2} \
      | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam

    samtools index ${sample_id}.sorted.bam
    """
}

process CALCULATE_COVERAGE {
    tag "$sample_id"
    publishDir "${params.outdir}/coverage", mode: 'copy'

    cpus 1
    memory "2 GB"
    conda "${projectDir}/envs/bwa_samtools.yml"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.coverage.txt")

    script:
    """
    samtools depth ${bam} > ${sample_id}.coverage.txt
    """
}

process PLOT_COVERAGE {
    tag "$sample_id"
    publishDir "${params.outdir}/coverage_plots", mode: 'copy'

    cpus 1
    memory "2 GB"
    conda "${projectDir}/envs/python_plot.yml"

    input:
    tuple val(sample_id), path(cov)

    output:
    tuple val(sample_id), path("${sample_id}.coverage.png")

    script:
    """
    python ${projectDir}/HW2/bin/plot_coverage.py \
      --input ${cov} \
      --output ${sample_id}.coverage.png \
      --sample ${sample_id}
    """
}

workflow QC_RAW {
    take:
    reads_ch

    main:
    raw_labeled_ch = reads_ch.map { sample_id, read1, read2 ->
        tuple(sample_id, "raw", read1, read2)
    }

    raw_qc_ch = FASTQC(raw_labeled_ch)

    emit:
    raw_qc_ch
}


workflow QC_TRIMMED {
    take:
    trimmed_reads_ch

    main:
    trimmed_labeled_ch = trimmed_reads_ch.map { sample_id, read1, read2 ->
        tuple(sample_id, "trimmed", read1, read2)
    }

    trimmed_qc_ch = FASTQC(trimmed_labeled_ch)

    emit:
    trimmed_qc_ch
}


workflow MAP_TO_REFERENCE {
    take:
    trimmed_reads_ch
    ref_index_ch

    main:
    bam_ch = MAP_READS(trimmed_reads_ch, ref_index_ch)

    emit:
    bam_ch
}

workflow COVERAGE_AND_PLOT {
    take:
    bam_ch

    main:
    cov_ch = CALCULATE_COVERAGE(bam_ch)
    plot_ch = PLOT_COVERAGE(cov_ch)

    emit:
    cov_ch
    plot_ch
}

workflow {
    reads_ch = Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .map { sample_id, reads -> tuple(sample_id, reads[0], reads[1]) }

    // 1. QC raw reads
    raw_qc_ch = QC_RAW(reads_ch)

    // 2. Trim reads
    trim_res = TRIM_READS(reads_ch)
    trimmed_reads_ch = trim_res.reads

    // 3. QC trimmed reads
    trimmed_qc_ch = QC_TRIMMED(trimmed_reads_ch)

    ref_index_ch = INDEX_REFERENCE(file(params.reference))
    bam_ch = MAP_TO_REFERENCE(trimmed_reads_ch, ref_index_ch)

    COVERAGE_AND_PLOT(bam_ch)
}