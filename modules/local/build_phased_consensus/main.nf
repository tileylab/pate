process BUILD_PHASED_CONSENSUS {
    tag "${meta.id}_${meta.locus}"
    label 'process_single'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(phase_out), path(vcf), path(reference_fasta), path(iupac_fasta)

    output:
    tuple val(meta), path("*.phased.fasta"), emit: phased_fasta
    tuple val(meta), path("*.genotype.fasta"), emit: genotype_fasta
    tuple val(meta), path("*.pickone.fasta"), emit: pickone_fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.locus}"
    def unique_flag = params.unique_only ? '--unique_only' : ''
    def dosage_flag = params.output_expected_dosage ? '--output_expected_dosage' : ''
    def genotype_mode = params.genotype_mode ?: 'consensus'
    """
    build_phased_consensus.py \\
        --phase_file ${phase_out} \\
        --vcf ${vcf} \\
        --reference_fasta ${reference_fasta} \\
        --iupac_fasta ${iupac_fasta} \\
        --sample_id ${meta.id} \\
        --locus ${meta.locus} \\
        --ploidy ${meta.ploidy} \\
        --genotype_mode ${genotype_mode} \\
        ${unique_flag} \\
        ${dosage_flag} \\
        --output_prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.locus}"
    """
    touch ${prefix}.phased.fasta
    touch ${prefix}.genotype.fasta
    touch ${prefix}.pickone.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
