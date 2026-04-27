process GENERATE_GENOME_BINS {
    tag "Bins: ${genome} (${bin_size}bp)"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/bedtools:2.31.1--7c4ce4cb07c09ee4' :
        'community.wave.seqera.io/library/bedtools:2.31.1--7c4ce4cb07c09ee4' }"

    input:
    path chrom_sizes    
    val genome          
    val bin_size        

    output:
    path "${out_filename}", emit: bins_file
    path "versions.yml"   , emit: versions

    script:
    out_filename = "${genome}_${bin_size}bp_bins.bed"
    """
    set -euo pipefail

    bedtools makewindows \\
        -g "${chrom_sizes}" \\
        -w "${bin_size}" \\
    | awk -v size="${bin_size}" 'BEGIN{OFS="\\t"} (\$3-\$2) == size' \\
    > "${out_filename}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    out_filename = "${genome}_${bin_size}bp_bins.bed"
    
    """
    touch "${out_filename}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
