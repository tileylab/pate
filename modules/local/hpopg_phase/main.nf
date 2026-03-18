process HPOPG_PHASE {
    tag "${meta.id}_${meta.locus}"
    label 'process_single'

    container 'gptiley/h-popg:latest'

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
    """
    java -jar /opt/H-PoPG.jar \\
        -b ${bam} \\
        -v ${vcf} \\
        -p ${meta.ploidy} \\
        -o ${prefix}.phase.out \\
        -d ${prefix}.phase.log \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.locus}"
    """
    touch ${prefix}.phase.out
    touch ${prefix}.phase.log
    """
}
