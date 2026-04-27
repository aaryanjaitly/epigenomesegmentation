process GET_GENOME_CHROM_SIZE {
    tag "Getting ${genome} chrom sizes"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://biocontainers/biocontainers:v1.2.0_cv1' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    val genome          

    output:
    path "${genome}.chrom.sizes", emit: genome_chrom_size_file
    path "versions.yml"         , emit: versions

    script:
    """
    # Execute external bash script
    get_chrom_sizes.sh "${genome}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: \$(curl --version | head -n 1 | awk '{print \$2}')
        awk: \$(awk --version | head -n 1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    touch "${genome}.chrom.sizes"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: \$(curl --version | head -n 1 | awk '{print \$2}')
        awk: \$(awk --version | head -n 1 | awk '{print \$3}')
    END_VERSIONS
    """
}