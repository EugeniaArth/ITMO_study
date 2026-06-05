process FILTER_VARIANTS {
    tag "${meta.id} (${meta.country})"
    label 'process_low'

    publishDir "${params.outdir}/filtered_variants", mode: 'copy'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.filtered.vcf.gz"), path("${meta.id}.filtered.vcf.gz.tbi")

    script:
    """
    bcftools view \\
      -f 'PASS' \\
      -i 'QUAL>=${params.min_qual}' \\
      ${vcf} \\
      -Oz -o ${meta.id}.filtered.vcf.gz

    tabix ${meta.id}.filtered.vcf.gz
    """

    stub:
    """
    echo '##fileformat=VCFv4.2' | gzip > ${meta.id}.filtered.vcf.gz
    touch ${meta.id}.filtered.vcf.gz.tbi
    """
}
