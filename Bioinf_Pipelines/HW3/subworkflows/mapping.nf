include { INDEX_REFERENCE } from '../modules/local/bwa/index/main'
include { MAP_READS } from '../modules/local/bwa/mem/main'

workflow MAP_TO_REFERENCE {
    take:
    trimmed_reads_ch
    reference

    main:
    ref_index_ch = INDEX_REFERENCE(reference)
    bam_ch = MAP_READS(trimmed_reads_ch, ref_index_ch)

    emit:
    bam_ch
    ref_index_ch
}
