process BAMTOOLS_SPLIT {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bamtools=2.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bamtools:2.5.2--hd03093a_3' :
        'quay.io/biocontainers/bamtools:2.5.2--hd03093a_3' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.REF_*.bam"), emit: bams
    tuple val("${task.process}"), val('bamtools'), eval('bamtools --version 2>&1 | head -1 | sed "s/.*: //"'), topic: versions, emit: versions_bamtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bamtools split \\
        -in ${bam} \\
        -reference \\
        ${args}
    """

    stub:
    """
    touch ${meta.id}.REF_locus1.bam
    """
}
