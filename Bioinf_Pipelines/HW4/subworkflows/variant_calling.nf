include { BCFTOOLS_MPILEUP } from '../modules/nf-core/bcftools/mpileup/main'

workflow VARIANT_CALLING {
    take:
    bam_ch        // tuple(meta, bam, bai)
    ref_index_ch

    main:
    ref_ch = ref_index_ch.map { ref, index_files ->
        def fai = index_files.find { it.name.endsWith('.fai') }
        tuple([id: 'reference'], ref, fai)
    }

    bam_for_vc = bam_ch.map { meta, bam, bai ->
        tuple(meta, bam, [], [])
    }

    vcf_ch = BCFTOOLS_MPILEUP(bam_for_vc, ref_ch, false)

    emit:
    vcf   = vcf_ch.vcf
    stats = vcf_ch.stats
}
