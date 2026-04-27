process MERGE_COUNTS {
    tag "Merge: ${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/bedtools:2.31.1--7c4ce4cb07c09ee4' :
        'community.wave.seqera.io/library/bedtools:2.31.1--7c4ce4cb07c09ee4' }"

    input:
    tuple val(meta), 
          val(meta_list_meth), 
          path(files_meth), 
          val(meta_list_hist), 
          path(files_hist) 
    val   genome
    path  ref            
    path  segment_coords 

    output:
    tuple val(meta), 
          path("${meta.id}_seg_tabs_counts.bed"), 
          emit: merged_counts
    path "versions.yml", emit: versions

    script:
    def sample_id = meta.id
    def wgbs_path = ""
    def nome_path = ""

    // Assign paths based on modality
    meta_list_meth.eachWithIndex { m, i ->
        def current_file = files_meth instanceof List ? files_meth[i].name : files_meth.name
        if (m.modality == 'WGBS') { wgbs_path = current_file } 
        else if (m.modality == 'NOMe-seq') { nome_path = current_file }
    }

    def histone_path = files_hist instanceof List ? files_hist[0] : files_hist

    """
    # Execute the external bash script and pass variables as arguments
    merge_counts.sh \\
        "${sample_id}" \\
        "${segment_coords}" \\
        "${wgbs_path}" \\
        "${nome_path}" \\
        "${histone_path}"

    # Generate version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        awk: \$(awk --version | head -n 1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    def sample_id = meta.id
    """
    touch "${sample_id}_seg_tabs_counts.bed"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        awk: \$(awk --version | head -n 1 | awk '{print \$3}')
    END_VERSIONS
    """
}