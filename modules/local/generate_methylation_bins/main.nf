process GENERATE_METHYLATION_BINS {

    tag "Generate methyl bins: ${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/bedtools:2.31.1--7c4ce4cb07c09ee4' :
        'community.wave.seqera.io/library/bedtools:2.31.1--7c4ce4cb07c09ee4' }"

    input:
    tuple val(meta),
          val(meth_modalities),
          path(meth_tabs)
    path  segment_coords

    output:
    tuple val(meta),
          path("${meta.id}_seg_tabs_counts.txt"),
          emit: merged_counts
    path "versions.yml", emit: versions

    script:
    def sample_id = meta.id

    def mods = meth_modalities instanceof List ? meth_modalities : [meth_modalities]
    def tabs = meth_tabs        instanceof List ? meth_tabs        : [meth_tabs]

    def wgbs_idx  = mods.findIndexOf { it == 'WGBS' }
    def wgbs_path = wgbs_idx >= 0 ? tabs[wgbs_idx] : null

    """
    set -euo pipefail

    echo "[GENERATE_METHYLATION_BINS] Processing: ${sample_id}"

    SORTED_SEG="${sample_id}_segments_sorted.bed"
    FINAL_OUT="${sample_id}_seg_tabs_counts.txt"

    bedtools sort -i "${segment_coords}" > "\$SORTED_SEG"

    if [[ -z "${wgbs_path ?: ''}" ]]; then
        echo "[GENERATE_METHYLATION_BINS] ERROR: WGBS modality required but not found"
        exit 1
    fi

    echo "[GENERATE_METHYLATION_BINS] Mapping WGBS: ${wgbs_path}"

    head -n 1 "${wgbs_path}" > header.txt
    tail -n +2 "${wgbs_path}" > wgbs_body.bed

    bedtools map \\
        -a "\$SORTED_SEG" \\
        -b wgbs_body.bed \\
        -c 4,5 \\
        -o median \\
        -null 0 |
    awk 'NR==FNR{print; next} {print}' header.txt - > "\$FINAL_OUT"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        awk: \$(awk --version | head -n 1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    def sample_id = meta.id

    """
    touch "${sample_id}_seg_tabs_counts.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        awk: \$(awk --version | head -n 1 | awk '{print \$3}')
    END_VERSIONS
    """
}