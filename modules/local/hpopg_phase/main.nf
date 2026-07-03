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
        # H-PoPG throws (e.g. IndexOutOfBoundsException) when a locus has variants but no
        # phasing-informative read coverage ("Total effective calls: 0"). Treat any H-PoPG
        # failure as an unphaseable locus: emit empty outputs and continue, mirroring the
        # no-variants skip above. Downstream build_phased_consensus.py handles empty phase
        # files by falling back to the reference/IUPAC sequence.
        set +e
        java -jar /app/H-PoPGv0.2.0.jar \\
            -b ${bam} \\
            -v \${LOCUS_VCF} \\
            -p ${meta.ploidy} \\
            -o ${prefix}.phase.out \\
            -d ${prefix}.phase.log ${args}
        HPOPG_EXIT=\$?
        set -e

        if [ "\$HPOPG_EXIT" -ne 0 ]; then
            echo "WARNING: H-PoPG exited \$HPOPG_EXIT for ${prefix} (likely no phasing-informative coverage). Emitting empty phase outputs." >&2
            # rm+touch (not a redirect) because the process shell keeps -C (noclobber) and
            # H-PoPG may have already created a partial phase.out before crashing.
            rm -f ${prefix}.phase.out ${prefix}.phase.log
            touch ${prefix}.phase.out ${prefix}.phase.log
        fi
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.locus}"
    """
    touch ${prefix}.phase.out
    touch ${prefix}.phase.log
    """
}
