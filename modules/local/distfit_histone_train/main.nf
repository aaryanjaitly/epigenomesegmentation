process DISTFIT_HISTONE_TRAIN {
    tag "${meta.id} | ${mark} | ${dist}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://aaryanjaitly/episegmix:new_plots' :
        'aaryanjaitly/episegmix:new_plots' }"

    beforeScript "export PATH=\$PATH:${projectDir}/bin/src:${projectDir}/bin/HMM/build; export PYTHONPATH=\$PYTHONPATH:/app/src:${projectDir}/bin/src"

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
    cat <<EOF > ${prefix}.yaml
states: 3
marker: 1
marker_spec:
  - name: ${mark}
    distribution: ${dist}
data: [\$(readlink -f ${histone_data})]
EOF

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

    # 4. Inject Topology directly
    python -c "
import json
with open('${prefix}.preinit.json') as f:
    hmm = json.load(f)
N = int(hmm['states'])
hmm['topology'] = {str(s+1): [s] for s in range(N)}
with open('${prefix}.init.json', 'w') as f:
    json.dump(hmm, f, indent=4)
"

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