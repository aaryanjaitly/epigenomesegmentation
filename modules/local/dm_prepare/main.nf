process DM_PREPARE {
    tag "Prep: ${meta.id} (State: ${state})"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://aaryanjaitly/episegmix:new_plots' :
        'aaryanjaitly/episegmix:new_plots' }"
    
    // Minimal paths: Mac bin for Conda, and PYTHONPATH so imports work everywhere
    beforeScript "export PATH=\$PATH:${projectDir}/bin/src; export PYTHONPATH=\$PYTHONPATH:/app/src:${projectDir}/bin/src"
    
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
    def dist_overrides = meta.distributions ? meta.distributions.collect { k, v -> "${k}:${v}" }.join(",") : ""

    """
    set -euo pipefail
    
    # [Internal setup logic remains same for YAML creation]
    HAS_METH=0
    if [ -s "${meth}" ]; then HAS_METH=1; fi

    HEADER=\$(head -n 1 "${histone}")
    FILE_MARKS=(\$(echo "\$HEADER" | cut -f4-))
    
    MARKER_SPEC=""
    for M in "\${FILE_MARKS[@]}"; do
        DIST="${dist_hist}"
        if [[ "${dist_overrides}" == *"\$M:"* ]]; then
             DIST=\$(echo "${dist_overrides}" | grep -o "\$M:[^,]*" | cut -d: -f2)
        fi
        MARKER_SPEC="\${MARKER_SPEC}  - name: \${M}\n    distribution: \${DIST}\n"
    done

    cat <<EOF > "${prefix}.yaml"
states: ${state}
marker: \${#FILE_MARKS[@]}
marker_spec:
\${MARKER_SPEC}data: [\$(pwd)/${histone}]
chr: ${chr_params}
EOF

    if [ "\$HAS_METH" -eq 1 ]; then
        METH_DIST="${dist_meth}"
        if [[ "${dist_overrides}" == *"WGBS:"* ]]; then
             METH_DIST=\$(echo "${dist_overrides}" | grep -o "WGBS:[^,]*" | cut -d: -f2)
        fi
        cat <<EOF >> "${prefix}.yaml"
dna_methylation: \${METH_DIST}
meth_data: [\$(pwd)/${meth}]
EOF
    fi

    get_counts.py \\
        -d "${prefix}.yaml" \\
        -c "${prefix}-train-counts.txt" \\
        -m "${prefix}-train-counts-meth.txt" \\
        -r "${prefix}-train-regions.txt"

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