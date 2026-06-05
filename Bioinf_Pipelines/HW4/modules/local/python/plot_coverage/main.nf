process PLOT_COVERAGE {
    tag "$sample_id"
    label 'process_low'

    publishDir "${params.outdir}/coverage_plots", mode: 'copy'

    input:
    tuple val(sample_id), path(cov)

    output:
    tuple val(sample_id), path("${sample_id}.coverage.png")

    script:
    """
    python ${projectDir}/bin/plot_coverage.py \\
      --input ${cov} \\
      --output ${sample_id}.coverage.png \\
      --sample ${sample_id}
    """
}
