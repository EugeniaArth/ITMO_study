/*
 * Variant calling subworkflow using nf-core BCFTOOLS_MPILEUP module
 * (bcftools mpileup + bcftools call)
 */
include { BCFTOOLS_MPILEUP } from '../modules/nf-core/bcftools/mpileup/main'

workflow VARIANT_CALLING {
    take:
    bam_ch
    ref_index_ch

    main:
    ref_ch = ref_index_ch
        .map { ref, index_files ->
            def fai = index_files.find { it.name.endsWith('.fai') }
            tuple([id: 'reference'], ref, fai)
        }
        .first()

    bam_for_vc = bam_ch.map { sample_id, bam, bai ->
        tuple([id: sample_id], bam, [], [])
    }

    vcf_ch = BCFTOOLS_MPILEUP(bam_for_vc, ref_ch, false)

    emit:
    vcf = vcf_ch.vcf
    stats = vcf_ch.stats
}
