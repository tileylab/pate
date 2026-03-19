process PICARD_FASTQTOSAM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.unmapped.bam"), emit: bam
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_picard

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_reads = meta.single_end ? "--FASTQ ${reads[0]}" : "--FASTQ ${reads[0]} --FASTQ2 ${reads[1]}"
    def avail_mem = 3072
    if (task.memory) {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        FastqToSam \\
        ${input_reads} \\
        --OUTPUT ${prefix}.unmapped.bam \\
        --SAMPLE_NAME ${meta.id} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.unmapped.bam
    """
}
