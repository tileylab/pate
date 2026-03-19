process HPOPG_PHASE {
    tag "${meta.id}_${meta.locus}"
    label 'process_single'

    container 'docker.io/gptiley/h-popg:latest'

    input:
    tuple val(meta), path(bam), path(bai), path(vcf)

    output:
    tuple val(meta), path("*.phase.out"), emit: phase_out
    tuple val(meta), path("*.phase.log"), emit: phase_log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.locus}"
    def vcf_unzipped = vcf.name.endsWith('.gz') ? vcf.name.replaceAll(/\.gz$/, '') : vcf.name
    """
    if [[ "${vcf}" == *.gz ]]; then
        gunzip -c ${vcf} > ${vcf_unzipped}
    else
        ln -s ${vcf} ${vcf_unzipped}
    fi

    # Filter VCF to only variants on the current locus
    # The whole-sample VCF contains variants from all loci, but the BAM is per-locus
    LOCUS_VCF="${prefix}.locus.vcf"
    grep '^#' ${vcf_unzipped} > \${LOCUS_VCF}
    grep -v '^#' ${vcf_unzipped} | awk -v locus="${meta.locus}" '\$1 == locus' >> \${LOCUS_VCF} || true

    # Count locus-specific variant records
    NVAR=\$(grep -cv '^#' \${LOCUS_VCF} || true)

    if [ "\$NVAR" -eq 0 ]; then
        echo "WARNING: No variants found for locus ${meta.locus}. Skipping H-PoPG phasing for ${prefix}." >&2
        touch ${prefix}.phase.out
        touch ${prefix}.phase.log
    else
        java -jar /app/H-PoPGv0.2.0.jar \\
            -b ${bam} \\
            -v \${LOCUS_VCF} \\
            -p ${meta.ploidy} \\
            -o ${prefix}.phase.out \\
            -d ${prefix}.phase.log \\
            ${args}
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.locus}"
    """
    touch ${prefix}.phase.out
    touch ${prefix}.phase.log
    """
}
