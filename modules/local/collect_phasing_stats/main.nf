process COLLECT_PHASING_STATS {
    tag "all_samples"
    label 'process_single'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    path(phased_fastas)
    path(vcfs)
    path(phase_outs)
    path(reference_fastas)
    val(sample_ploidies)

    output:
    path("*.phasingSummary.txt"), emit: per_sample_stats
    path("averagePhasingStats.txt"), emit: avg_stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ploidy_json = groovy.json.JsonOutput.toJson(sample_ploidies)
    """
    collect_phasing_stats.py \\
        --phased_dir . \\
        --vcf_dir . \\
        --phase_dir . \\
        --reference_dir . \\
        --ploidy_json '${ploidy_json}' \\
        --output_dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch sample1.phasingSummary.txt
    touch averagePhasingStats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
