process DNA_REPORT {
    tag "Report: ${meta.id}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://aaryanjaitly/episegmix:new_plots' :
        'aaryanjaitly/episegmix:new_plots' }"

    beforeScript """
        export PATH=\$PATH:${projectDir}/bin/workflowDNAMethylation:${projectDir}/bin/src
        export PYTHONPATH=\$PYTHONPATH:/app/src:${projectDir}/bin/src
    """

    input:
    tuple val(meta), path(config), path(model), path(bed), path(seg_txt)

    output:
    path "plots/${meta.id}/*.png",  emit: plots
    path "plots/${meta.id}/*.html", emit: html_report
    path "versions.yml",            emit: versions

    script:
    """
    set -euo pipefail
    export MPLCONFIGDIR=\$(pwd)
    mkdir -p plots/${meta.id}
    
    # 1. GENERATE RESULTS PLOTS
    results.py \\
        -c "${seg_txt}" \\
        -j "${model}" \\
        -e "plots/${meta.id}/${meta.id}-meanEmission.png" \\
        -n "plots/${meta.id}/${meta.id}-normEmission.png" \\
        -t "plots/${meta.id}/${meta.id}-transitionMatrix.png" \\
        -m "plots/${meta.id}/${meta.id}-stateMembership.png" \\
        -l "plots/${meta.id}/${meta.id}-stateLength.png" \\
        -s viterbi \\
        -d ${bed}

    # 2. GENERATE STATE COLORS
    plot_state_colors.py \\
        -d ${bed} \\
        -o "plots/${meta.id}/${meta.id}-stateColors.png"

    # 3. RUN DNA-SPECIFIC REPORT
    segmentation_report_dna.sh \\
        -n "${meta.id}" \\
        -o "plots/${meta.id}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p plots/${meta.id}
    touch plots/${meta.id}/dummy.png
    touch plots/${meta.id}/dummy.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """
}