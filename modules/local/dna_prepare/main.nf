process DNA_PREPARE {
    tag "Prep: ${meta.id} (State: ${state})"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://aaryanjaitly/episegmix:new_plots' :
        'aaryanjaitly/episegmix:new_plots' }"

    beforeScript """
        export PATH=\$PATH:${projectDir}/bin/src; 
        export PYTHONPATH=\$PYTHONPATH:/app/src:${projectDir}/bin/src
    """

    input:
    tuple val(meta), path(meth), val(state)

    output:
    tuple val(meta), path("${meta.id}.yaml"),             emit: config
    tuple val(meta), path("${meta.id}.trainCounts.txt"),  emit: train_counts
    tuple val(meta), path("${meta.id}.trainregions.txt"), emit: train_regions
    tuple val(meta), path("${meta.id}.counts.txt"),       emit: counts
    tuple val(meta), path("${meta.id}.regions.txt"),      emit: regions
    path "versions.yml",                                  emit: versions

    script:
    def prefix = meta.id
    def meth_files = (meth instanceof List ? meth : [meth]).join(' ')
    def dist_meth = params.dist_methyl.toString()
    def dist_overrides = meta.distributions ? meta.distributions.collect { k, v -> "${k}:${v}" }.join(",") : ""

    """
    set -euo pipefail
    
    # 1. Create YAML config via external script
    create_dna_yaml.py \\
        --prefix "${prefix}" \\
        --meth ${meth_files} \\
        --state "${state}" \\
        --dist_meth "${dist_meth}" \\
        --overrides "${dist_overrides}"

    # 2. Generate Counts
    get_meth_counts.py \\
        -d "\$(pwd)/${prefix}.yaml" \\
        -c "\$(pwd)/${prefix}.trainCounts.txt" \\
        -r "\$(pwd)/${prefix}.trainregions.txt" \\
        -C "\$(pwd)/${prefix}.counts.txt" \\
        -R "\$(pwd)/${prefix}.regions.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.yaml"
    touch "${prefix}.trainCounts.txt"
    touch "${prefix}.trainregions.txt"
    touch "${prefix}.counts.txt"
    touch "${prefix}.regions.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
END_VERSIONS
    """
}