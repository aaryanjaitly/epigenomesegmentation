process DM_PREPARE {
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
    tuple val(meta), path(histone), path(meth), val(state)
    val chr_params

    output:
    tuple val(meta), path("${meta.id}.yaml"),                  emit: config
    tuple val(meta), path("${meta.id}-train-counts.txt"),      emit: train_counts
    tuple val(meta), path("${meta.id}-train-regions.txt"),     emit: regions
    tuple val(meta), path("${meta.id}-train-counts-meth.txt"), emit: train_meth
    path "versions.yml",                                       emit: versions

    script:
    def prefix    = "${meta.id}"
    def dist_hist = params.dist_histone.toString()
    def dist_meth = params.dist_methyl.toString()
    
    // Safely format the overrides map into a comma-separated string for the python script
    def dist_overrides = meta.distributions ? meta.distributions.collect { k, v -> "${k}:${v}" }.join(",") : ""

    """
    set -euo pipefail
    
    # 1. Generate the YAML config via external script
    create_dm_yaml.py \\
        --prefix "${prefix}" \\
        --histone "${histone}" \\
        --meth "${meth}" \\
        --state "${state}" \\
        --chr_params "${chr_params}" \\
        --dist_hist "${dist_hist}" \\
        --dist_meth "${dist_meth}" \\
        --overrides "${dist_overrides}"

    # 2. Get Counts using the newly generated YAML
    get_counts.py \\
        -d "${prefix}.yaml" \\
        -c "${prefix}-train-counts.txt" \\
        -m "${prefix}-train-counts-meth.txt" \\
        -r "${prefix}-train-regions.txt"

    # 3. Capture versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    touch "${prefix}.yaml"
    touch "${prefix}-train-counts.txt"
    touch "${prefix}-train-regions.txt"
    touch "${prefix}-train-counts-meth.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
END_VERSIONS
    """
}