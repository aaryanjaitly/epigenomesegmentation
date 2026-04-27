process DISTFIT_HISTONE_ASSESS {
    tag "${meta.id} - Evaluate"
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
    tuple val(meta), path(histone_data), path(models)
    path original_csv 

    output:
    tuple val(meta), path("DISTFIT_${meta.id}_samplesheet.csv"), emit: updated_samplesheet
    path "versions.yml", emit: versions 

    script:
    def model_files = models.join(' ')
    """
    set -euo pipefail
    export MPLCONFIGDIR=\$(pwd)

    # 1. Compile LogLikelihood if missing
    if ! command -v LogLikelihood &> /dev/null || ! LogLikelihood 2>&1 | grep -qi "usage"; then
        bash "${projectDir}/bin/build.sh"
        export PATH="\$PATH:${projectDir}/bin/HMM/build"
    fi

    # 2. Extract markers
    MARKERS=\$(head -n 1 ${histone_data} | awk '{for(i=4;i<=NF;++i) print \$i}')

    # 3. Generate counts per marker
    for MARK in \$MARKERS; do
        mkdir -p temp_counts/\$MARK
        echo -e "states: 3\\nmarker: 1\\nmarker_spec:\\n  - name: \$MARK\\n    distribution: NBI\\ndata: [\$(readlink -f ${histone_data})]" > dummy_\${MARK}.yaml
        get_counts_for_all.py -d dummy_\${MARK}.yaml -o temp_counts/\$MARK
    done

    # 4. Compute log likelihoods
    touch log_likelihoods.txt
    for model in ${model_files}; do
        fname=\$(basename "\$model" .final-model.json)
        dist=\$(echo "\$fname" | cut -d'-' -f1)
        mark=\$(echo "\$fname" | awk -F'-model-' '{print \$2}')

        for cfile in temp_counts/\${mark}/counts*.txt; do
            [ -f "\$cfile" ] || continue
            rfile="temp_counts/\${mark}/regions_\${cfile##*counts_}"
            
            if [ -f "\$rfile" ]; then
                score=\$(LogLikelihood -m "\$model" -c "\$cfile" -r "\$rfile" || echo "NaN")
                echo -e "\$dist\\t\$mark\\t\$score" >> log_likelihoods.txt
            fi
        done
    done

    # 5. Update samplesheet via external script
    update_distfit_samplesheet.py \\
        -l log_likelihoods.txt \\
        -c ${original_csv} \\
        -s ${meta.id} \\
        -o DISTFIT_${meta.id}_samplesheet.csv

    # 6. Capture versions 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
END_VERSIONS
    """

    stub:
    """
    touch "DISTFIT_${meta.id}_samplesheet.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
END_VERSIONS
    """
}