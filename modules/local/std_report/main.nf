process STD_REPORT {
    tag "Report: ${meta.id}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://aaryanjaitly/episegmix:new_plots' :
        'aaryanjaitly/episegmix:new_plots' }"

    // Environmental setup for both Conda (local bin) and Docker (system PATH)
    beforeScript """
        export PATH=\$PATH:${projectDir}/bin/workflow:${projectDir}/bin/src
        export PYTHONPATH=\$PYTHONPATH:/app/src:${projectDir}/bin/src
    """

    input:
    tuple val(meta), path(config), path(model), path(bed), path(seg_txt)

    output:
    path "plots/${meta.id}/*.png" , emit: plots
    path "plots/${meta.id}/*.html", emit: html_report
    path "versions.yml"           , emit: versions

    script:
    """
    set -euo pipefail
    export MPLCONFIGDIR=\$(pwd)
    
    mkdir -p plots/${meta.id}

    # 1. GENERATE DATA STATISTICS PLOTS
    # Found in bin/src (Conda) or /usr/local/bin (Docker)
    plot_statistics.py \\
        -d "${config}" \\
        -p "plots/${meta.id}/${meta.id}-histogram.png" \\
        -c "plots/${meta.id}/${meta.id}-correlation.png" \\
        -m "plots/${meta.id}/${meta.id}-methylation-density.png"

    # 2. GENERATE MODEL EMISSION & TRANSITION PLOTS
    results.py \\
        -c "${seg_txt}" \\
        -j "${model}" \\
        -e "plots/${meta.id}/${meta.id}-meanEmission-viterbi.png" \\
        -t "plots/${meta.id}/${meta.id}-transitionMatrix.png" \\
        -m "plots/${meta.id}/${meta.id}-stateMembership-viterbi.png" \\
        -l "plots/${meta.id}/${meta.id}-stateLength-viterbi.png" \\
        -s viterbi \\
        -n "plots/${meta.id}/${meta.id}-normEmission-viterbi.png" \\
        -d "${bed}"

    # 3. GENERATE STATE DISTRIBUTION HISTOGRAMS
    plot_state_histograms.py \\
        -c "${seg_txt}" \\
        -j "${model}" \\
        -a "plots/${meta.id}/${meta.id}-stateDistribution-viterbi.png" \\
        -s viterbi

    # 4. GENERATE GENOMIC STATE COLOR MAPS
    plot_state_colors.py \\
        -d "${bed}" \\
        -o "plots/${meta.id}/${meta.id}-state-colors.png"

    # 5. ASSEMBLE FINAL HTML REPORT (Standard Mode)
    # Found in bin/workflow (Conda) or /usr/local/bin (Docker)
    segmentation_report_std.sh \\
        -n "${meta.id}" \\
        -o "plots/${meta.id}" \\
        -i true \\
        -c viterbi

    # 6. TRACK VERSIONING
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