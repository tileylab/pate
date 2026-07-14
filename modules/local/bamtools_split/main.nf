process BAMTOOLS_SPLIT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d7/d7e24dc1e4d93ca4d3a76a78d4c834a7be3985b0e1e56fddd61662e047863a8a/data' :
        'community.wave.seqera.io/library/bwa_htslib_samtools:83b50ff84ead50d0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.REF_*.bam"), path("*.REF_*.bam.bai"), emit: bams
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_bamtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: bam.baseName
    """
    # Index the input BAM if no index exists
    if [ ! -f ${bam}.bai ]; then
        samtools index ${bam}
    fi

    # Split BAM by reference contig (equivalent to bamtools split -reference) and index each
    # per-locus BAM in the same loop (H-PoPG needs the .bai for random access). Indexing here
    # avoids a separate per-locus index job fan-out.
    for ref in \$(samtools idxstats ${bam} | cut -f1 | grep -v '^\\*'); do
        samtools view -b ${bam} "\${ref}" >| "${prefix}.REF_\${ref}.bam"
        samtools index "${prefix}.REF_\${ref}.bam"
    done
    """

    stub:
    """
    touch ${meta.id}.REF_locus1.bam
    touch ${meta.id}.REF_locus1.bam.bai
    """
}
