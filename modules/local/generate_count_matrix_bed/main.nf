process GENERATE_COUNT_MATRIX_BED {
    tag "Counts BED: ${meta.id} | ${meta.modality}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://biocontainers/biocontainers:v1.2.0_cv1' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    tuple val(meta), path(file_path)

    output:
    tuple val(meta), path("${meta.id}_${meta.modality}/${meta.id}_${meta.modality}.tab"), emit: counts_bed
    path "versions.yml",                                                                  emit: versions

    script:
    def sample_id = meta.id
    def modality  = meta.modality
    
    """
    # Execute external bash script
    generate_count_matrix_bed.sh "${sample_id}" "${modality}" "${file_path}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n 1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    def sample_id = meta.id
    def modality  = meta.modality

    """
    mkdir -p "${sample_id}_${modality}"
    touch "${sample_id}_${modality}/${sample_id}_${modality}.tab"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n 1 | awk '{print \$3}')
    END_VERSIONS
    """
}