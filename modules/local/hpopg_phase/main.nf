process HPOPG_PHASE {
    tag "${meta.id}"
    label 'process_single'

    container 'docker.io/gptiley/h-popg:latest'

    input:
    tuple val(meta), path(bams), path(bais), path(vcf)

    output:
    tuple val(meta), path("*.phase.out"), emit: phase_out
    tuple val(meta), path("*.phase.log"), emit: phase_log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def vcf_unzipped = vcf.name.endsWith('.gz') ? vcf.name.replaceAll(/\.gz$/, '') : vcf.name
    """
    # Decompress the per-sample filtered VCF once (idempotent under -C noclobber)
    if [[ "${vcf}" == *.gz ]]; then
        gunzip -c ${vcf} >| ${vcf_unzipped}
    else
        ln -sf ${vcf} ${vcf_unzipped}
    fi

    # Phase every per-locus BAM produced by BAMTOOLS_SPLIT for this sample.
    # One task per sample loops over all loci to avoid a per-locus job fan-out.
    for bam in *.REF_*.bam; do
        locus=\${bam##*.REF_}
        locus=\${locus%.bam}
        prefix="${meta.id}.\${locus}"
        LOCUS_VCF="\${prefix}.locus.vcf"

        # Filter the whole-sample VCF to only variants on this locus
        grep '^#' ${vcf_unzipped} >| \${LOCUS_VCF}
        grep -v '^#' ${vcf_unzipped} | awk -v locus="\${locus}" '\$1 == locus' >> \${LOCUS_VCF} || true

        NVAR=\$(grep -cv '^#' \${LOCUS_VCF} || true)

        if [ "\$NVAR" -eq 0 ]; then
            echo "WARNING: No variants found for locus \${locus}. Skipping H-PoPG phasing for \${prefix}." >&2
            touch \${prefix}.phase.out \${prefix}.phase.log
            continue
        fi

        # H-PoPG throws (e.g. IndexOutOfBoundsException) when a locus has variants but no
        # phasing-informative read coverage ("Total effective calls: 0"). Treat any H-PoPG
        # failure as an unphaseable locus: emit empty outputs and continue, mirroring the
        # no-variants skip above. Downstream build_phased_consensus.py handles empty phase
        # files by falling back to the reference/IUPAC sequence.
        set +e
        java -jar /app/H-PoPGv0.2.0.jar \\
            -b \${bam} \\
            -v \${LOCUS_VCF} \\
            -p ${meta.ploidy} \\
            -o \${prefix}.phase.out \\
            -d \${prefix}.phase.log ${args}
        HPOPG_EXIT=\$?
        set -e

        if [ "\$HPOPG_EXIT" -ne 0 ]; then
            echo "WARNING: H-PoPG exited \$HPOPG_EXIT for \${prefix} (likely no phasing-informative coverage). Emitting empty phase outputs." >&2
            # rm+touch (not a redirect) because the process shell keeps -C (noclobber) and
            # H-PoPG may have already created a partial phase.out before crashing.
            rm -f \${prefix}.phase.out \${prefix}.phase.log
            touch \${prefix}.phase.out \${prefix}.phase.log
        fi
    done
    """

    stub:
    """
    touch ${meta.id}.stub.phase.out
    touch ${meta.id}.stub.phase.log
    """
}
