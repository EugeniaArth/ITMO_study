#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { QC_AND_TRIM } from './subworkflows/qc_and_trim'
include { MAP_TO_REFERENCE } from './subworkflows/mapping'
include { COVERAGE_AND_PLOT } from './subworkflows/coverage'
include { VARIANT_CALLING } from './subworkflows/variant_calling'

workflow {
    reads_ch = Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .map { sample_id, reads -> tuple(sample_id, reads[0], reads[1]) }

    qc_trim = QC_AND_TRIM(reads_ch)

    mapping = MAP_TO_REFERENCE(qc_trim.trimmed_reads_ch, file(params.reference))

    COVERAGE_AND_PLOT(mapping.bam_ch)

    VARIANT_CALLING(mapping.bam_ch, mapping.ref_index_ch)
}
