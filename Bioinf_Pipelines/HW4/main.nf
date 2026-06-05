#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { QC_AND_TRIM } from './subworkflows/qc_and_trim'
include { MAP_TO_REFERENCE } from './subworkflows/mapping'
include { COVERAGE_AND_PLOT } from './subworkflows/coverage'
include { VARIANT_CALLING as VARIANT_CALLING_USA    } from './subworkflows/variant_calling'
include { VARIANT_CALLING as VARIANT_CALLING_FRANCE } from './subworkflows/variant_calling'
include { FILTER_VARIANTS } from './modules/local/bcftools/filter/main'

workflow {
    
    // 1. Read all samples from CSV into ONE channel (map)
 
    samples_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample, country: row.country]
            tuple(meta, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true))
        }

    reads_ch = samples_ch.map { meta, r1, r2 -> tuple(meta.id, r1, r2) }

    // common steps for all samples (run once)
    qc_trim  = QC_AND_TRIM(reads_ch)
    mapping  = MAP_TO_REFERENCE(qc_trim.trimmed_reads_ch, file(params.reference, checkIfExists: true))
    COVERAGE_AND_PLOT(mapping.bam_ch)

    // attach metadata back to BAM (combine)
    meta_ch = samples_ch.map { meta, r1, r2 -> tuple(meta.id, meta) }
    bam_with_meta = mapping.bam_ch
        .combine(meta_ch, by: 0)
        .map { sample_id, bam, bai, meta -> tuple(meta, bam, bai) }

    
    // 2. Split channel by country field (branch)
    
    by_country = bam_with_meta.branch { meta, bam, bai ->
        usa:    meta.country == 'USA'
        france: meta.country == 'France'
    }

    
    // 3. process each group — variant calling
   
    usa_vcfs    = VARIANT_CALLING_USA(by_country.usa, mapping.ref_index_ch)
    france_vcfs = VARIANT_CALLING_FRANCE(by_country.france, mapping.ref_index_ch)

    // 4. Join groups back into ONE channel (mix)

    all_vcfs = usa_vcfs.vcf.mix(france_vcfs.vcf)


    // 5. Final analysis — filter variants (stub for testing)

    FILTER_VARIANTS(all_vcfs)
}
