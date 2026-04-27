process DISTFIT_HISTONE_ASSESS {
    tag "${meta.id} - Evaluate"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://aaryanjaitly/episegmix:new_plots' :
        'aaryanjaitly/episegmix:new_plots' }"

    beforeScript "export PATH=\$PATH:${projectDir}/bin/src:${projectDir}/bin/HMM/build; export PYTHONPATH=\$PYTHONPATH:/app/src:${projectDir}/bin/src"

    input:
    tuple val(meta), path(histone_data), path(models)
    path original_csv 

    output:
    tuple val(meta), path("DISTFIT_${meta.id}_samplesheet.csv"), emit: updated_samplesheet
    path "versions.yml", emit: versions 

    script:
    """
    set -euo pipefail
    export MPLCONFIGDIR=\$(pwd)

    # Auto-build for Conda: Checks if LogLikelihood exists AND can run successfully
    if ! command -v LogLikelihood &> /dev/null || ! LogLikelihood 2>&1 | grep -qi "usage"; then
        echo "LogLikelihood not found or incompatible. Running build.sh..."
        bash "${projectDir}/bin/build.sh"
        export PATH="\$PATH:${projectDir}/bin/HMM/build"
    fi

    RAW_HEADER=\$(head -n 1 ${histone_data})
    export HEADER_STR="\$RAW_HEADER"
    export HISTONE_ABS=\$(readlink -f ${histone_data})

    MARKERS=\$(python -c "import os; print(' '.join(os.environ.get('HEADER_STR', '').strip().split()[3:]))")

    mkdir -p temp_counts
    touch log_likelihoods.txt

    # 1. Generate counts per marker
    for MARK in \$MARKERS; do
        cat <<EOF > dummy_\${MARK}.yaml
states: 3
marker: 1
marker_spec:
  - name: \${MARK}
    distribution: NBI
data: [\$HISTONE_ABS]
EOF
        
        mkdir -p temp_counts/\${MARK}
        get_counts_for_all.py -d dummy_\${MARK}.yaml -o temp_counts/\${MARK}
    done

    # 2. Compute log likelihoods
    for model in ${models.join(' ')}; do
        fname=\$(basename "\$model" .final-model.json)
        dist=\${fname%%-*}
        mark=\${fname#*-model-}

        for cfile in temp_counts/\${mark}/counts*.txt; do
            [ -f "\$cfile" ] || continue
            base=\$(basename \$cfile)
            suffix="\${base#counts_}"
            rfile="temp_counts/\${mark}/regions_\$suffix"

            if [ -f "\$rfile" ]; then
                echo -n -e "\$dist\t\$mark\t" >> log_likelihoods.txt
                LogLikelihood \\
                    -m "\$model" \\
                    -c "\$cfile" \\
                    -r "\$rfile" >> log_likelihoods.txt || true
                echo "" >> log_likelihoods.txt
            fi
        done
    done

    # 3. Create final Samplesheet using Pandas
    python - <<EOF
import pandas as pd
import os

best_dists = {}
if os.path.exists('log_likelihoods.txt') and os.path.getsize('log_likelihoods.txt') > 0:
    data = []
    with open('log_likelihoods.txt') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                try:
                    score = float(parts[-1])
                    data.append({'Mark': parts[1], 'Dist': parts[0], 'Score': score})
                except ValueError:
                    continue
    if data:
        df = pd.DataFrame(data)
        best_models = df.loc[df.groupby('Mark')['Score'].idxmax()]
        for _, row in best_models.iterrows():
            best_dists[row['Mark']] = row['Dist']

orig_df = pd.read_csv('${original_csv}')
sample_df = orig_df[orig_df['sample_id'].astype(str) == "${meta.id}"].copy()

for idx, row in sample_df.iterrows():
    mark = row['epigenetic_mark']
    if mark in best_dists:
        sample_df.at[idx, 'distribution'] = best_dists[mark]

sample_df.to_csv("DISTFIT_${meta.id}_samplesheet.csv", index=False)
EOF

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