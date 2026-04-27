process STD_DECODE {
    tag "Decode: ${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://aaryanjaitly/episegmix:new_plots' :
        'aaryanjaitly/episegmix:new_plots' }"

    beforeScript "export PATH=\$PATH:${projectDir}/bin/src:${projectDir}/bin/HMM/build; export PYTHONPATH=\$PYTHONPATH:/app/src:${projectDir}/bin/src"

    input:
    tuple val(meta), path(config), path(model)

    output:
    tuple val(meta), path("segmentation/${meta.id}.bed.gz"), emit: bed
    tuple val(meta), path("segmentation/${meta.id}.txt")   , emit: seg_txt
    tuple val(meta), path("states_${meta.id}")             , emit: states_dir
    tuple val(meta), path("counts_${meta.id}")             , emit: counts_dir
    path "versions.yml"                                    , emit: versions

    script:
    """
    set -euo pipefail
    
    if ! command -v HMMChromSeg &> /dev/null; then
        bash "${projectDir}/bin/build.sh"
    fi
    
    mkdir -p counts_${meta.id} states_${meta.id} segmentation

    get_counts_for_all.py -d "${config}" -o "counts_${meta.id}"

    for file in \$(find "counts_${meta.id}" -name "counts*" -type f); do
        name=\${file##*/counts_}
        HMMChromSeg \\
            -m "${model}" \\
            -v "states_${meta.id}/viterbi_\$name" \\
            -d "states_${meta.id}/posterior_\$name" \\
            -c "\$file" \\
            -x "counts_${meta.id}/meth_\$name" \\
            -r "counts_${meta.id}/regions_\$name" \\
            -p ${task.cpus}
    done

    segmentation_to_bed.py -d "${config}" -i "states_${meta.id}/" -o "segmentation/" -c viterbi

    mv segmentation/*.bed.gz segmentation/${meta.id}.bed.gz
    mv segmentation/*.txt    segmentation/${meta.id}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p segmentation states_${meta.id} counts_${meta.id}
    touch segmentation/${meta.id}.bed.gz
    touch segmentation/${meta.id}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """
}