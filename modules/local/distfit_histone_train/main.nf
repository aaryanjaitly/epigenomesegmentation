process DISTFIT_HISTONE_TRAIN {
    tag "${meta.id} | ${mark} | ${dist}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://aaryanjaitly/episegmix:new_plots' :
        'aaryanjaitly/episegmix:new_plots' }"

    beforeScript """
        export PATH=\$PATH:${projectDir}/bin/src:${projectDir}/bin/HMM/build; 
        export PYTHONPATH=\$PYTHONPATH:/app/src:${projectDir}/bin/src
    """

    input:
    tuple val(meta), path(histone_data), val(mark), val(dist)

    output:
    tuple val(meta), path("*final-model.json"), emit: model
    tuple val(meta), val(mark), val(dist),      emit: meta_info 
    path "versions.yml",                        emit: versions

    script:
    def prefix = "${dist}-model-${mark}"

    """
    set -euo pipefail
    export MPLCONFIGDIR=\$(pwd)

    if ! command -v TopologyHMM &> /dev/null; then
        bash "${projectDir}/bin/build.sh"
    fi

    # 1. Config (Strictly Histone)
    echo -e "states: 3\\nmarker: 1\\nmarker_spec:\\n  - name: ${mark}\\n    distribution: ${dist}\\ndata: [\$(readlink -f ${histone_data})]" > ${prefix}.yaml

    # 2. Get Counts
    get_counts.py \\
        -d "${prefix}.yaml" \\
        -m "${prefix}.trainCountsMeth.txt" \\
        -c "${prefix}.trainCounts.txt" \\
        -r "${prefix}.trainRegions.txt"

    # 3. Init HMM
    init_HMM.py \\
        -d "${prefix}.trainCounts.txt" \\
        -e "${prefix}.trainCountsMeth.txt" \\
        -m "${prefix}.yaml" \\
        -j "${prefix}.preinit.json"

    # 4. Inject Topology via external script
    inject_topology.py \\
        -i "${prefix}.preinit.json" \\
        -o "${prefix}.init.json"

    # 5. Train HMM
    TopologyHMM \\
        -t -n ${params.adjustment} \\
        -m "${prefix}.init.json" \\
        -o "${prefix}.final-model.json" \\
        -c "${prefix}.trainCounts.txt" \\
        -r "${prefix}.trainRegions.txt" \\
        -i ${params.iter} -p ${task.cpus} -e ${params.epsilon}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
END_VERSIONS
    """

    stub:
    def prefix = "${dist}-model-${mark}"

    """
    touch "${prefix}.final-model.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
END_VERSIONS
    """
}