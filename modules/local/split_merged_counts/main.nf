process SPLIT_MERGED_COUNTS {
    tag "Split: ${meta.id}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://biocontainers/biocontainers:v1.2.0_cv1' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    tuple val(meta), path(merged_bed)

    output:
    tuple val(meta), 
          path("${meta.id}_Histone_Input.txt"), 
          path("${meta.id}_Meth_Input.txt"), 
          emit: inputs
    path "versions.yml", emit: versions

    script:
    def sample_id = meta.id
    """
    # Execute external bash script
    split_merge_counts.sh "${sample_id}" "${merged_bed}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n 1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    def sample_id = meta.id
    """
    touch "${sample_id}_Histone_Input.txt"
    touch "${sample_id}_Meth_Input.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n 1 | awk '{print \$3}')
    END_VERSIONS
    """
}