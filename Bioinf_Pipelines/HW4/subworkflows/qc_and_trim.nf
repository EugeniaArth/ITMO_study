include { FASTQC } from '../modules/local/fastqc/main'
include { TRIM_READS } from '../modules/local/fastp/main'

workflow QC_RAW {
    take:
    reads_ch

    main:
    raw_labeled_ch = reads_ch.map { sample_id, read1, read2 ->
        tuple(sample_id, 'raw', read1, read2)
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
        tuple(sample_id, 'trimmed', read1, read2)
    }
    trimmed_qc_ch = FASTQC(trimmed_labeled_ch)

    emit:
    trimmed_qc_ch
}

workflow QC_AND_TRIM {
    take:
    reads_ch

    main:
    raw_qc_ch = QC_RAW(reads_ch)

    trim_res = TRIM_READS(reads_ch)
    trimmed_reads_ch = trim_res.reads

    trimmed_qc_ch = QC_TRIMMED(trimmed_reads_ch)

    emit:
    raw_qc_ch
    trimmed_reads_ch
    trimmed_qc_ch
}
