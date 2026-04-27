process STD_TRAIN {
    tag "Train: ${meta.id}"
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://aaryanjaitly/episegmix:new_plots' :
        'aaryanjaitly/episegmix:new_plots' }"

    beforeScript """
        export PATH=\$PATH:${projectDir}/bin/src:${projectDir}/bin/HMM/build; 
        export PYTHONPATH=\$PYTHONPATH:/app/src:${projectDir}/bin/src
    """

    input:
    tuple val(meta), path(config), path(counts), path(regions), path(meth_counts)

    output:
    tuple val(meta), path("final-model-${meta.id}.json"), emit: model
    path "${meta.id}.log"                               , emit: log
    path "versions.yml"                                 , emit: versions

    script:
    def meth_arg = meth_counts ? "-e ${meth_counts}" : ""
    def train_meth_arg = meth_counts ? "-x ${meth_counts}" : ""

    """
    set -euo pipefail
    
    # Auto-build for Conda: Checks if HMMChromSeg exists AND can run successfully
    if ! command -v HMMChromSeg &> /dev/null || ! HMMChromSeg 2>&1 | grep -q "Usage"; then
        bash "${projectDir}/bin/build.sh"
    fi
    
    init_HMM.py \\
        -d "${counts}" \\
        ${meth_arg} \\
        -m "${config}" \\
        -j "${meta.id}-init.json"
    
    HMMChromSeg \\
        -t \\
        -m "${meta.id}-init.json" \\
        -o "final-model-${meta.id}.json" \\
        -c "${counts}" \\
        ${train_meth_arg} \\
        -r "${regions}" \\
        -i ${params.iter} -e ${params.epsilon} \\
        -p ${task.cpus} \\
        &> "${meta.id}.log"

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | awk '{print \$2}')
END_VERSIONS
    """

    stub:
    """
    touch "final-model-${meta.id}.json"
    touch "${meta.id}.log"

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | awk '{print \$2}')
END_VERSIONS
    """
}