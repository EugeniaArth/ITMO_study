include { CALCULATE_COVERAGE } from '../modules/local/samtools/depth/main'
include { PLOT_COVERAGE } from '../modules/local/python/plot_coverage/main'

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
