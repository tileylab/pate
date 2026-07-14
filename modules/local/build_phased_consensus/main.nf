process BUILD_PHASED_CONSENSUS {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(phase_outs), path(vcf), path(reference_fasta), path(iupac_fasta)

    output:
    tuple val(meta), path("*.phased.fasta"), emit: phased_fasta
    tuple val(meta), path("*.genotype.fasta"), emit: genotype_fasta
    tuple val(meta), path("*.pickone.fasta"), emit: pickone_fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def unique_flag = params.unique_only ? '--unique_only' : ''
    def dosage_flag = params.output_expected_dosage ? '--output_expected_dosage' : ''
    def genotype_mode = params.genotype_mode ?: 'consensus'
    """
    # One task per sample: build consensus for every per-locus phase output.
    # Locus is recovered from the phase-file name (\${meta.id}.\${locus}.phase.out).
    for pf in *.phase.out; do
        locus=\${pf#${meta.id}.}
        locus=\${locus%.phase.out}

        build_phased_consensus.py \\
            --phase_file \${pf} \\
            --vcf ${vcf} \\
            --reference_fasta ${reference_fasta} \\
            --iupac_fasta ${iupac_fasta} \\
            --sample_id ${meta.id} \\
            --locus \${locus} \\
            --ploidy ${meta.ploidy} \\
            --genotype_mode ${genotype_mode} \\
            ${unique_flag} \\
            ${dosage_flag} \\
            --output_prefix ${meta.id}.\${locus}
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.stub.phased.fasta
    touch ${meta.id}.stub.genotype.fasta
    touch ${meta.id}.stub.pickone.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
