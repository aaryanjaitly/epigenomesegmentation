process CLEAN_AND_INDEX_BAM {
    tag "Clean & Index: ${meta.id} | ${meta.modality}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/samtools:1.22.1--eccb42ff8fb55509' :
        'community.wave.seqera.io/library/samtools:1.22.1--eccb42ff8fb55509' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), 
          path("*.nochr.bam"), 
          path("*.nochr.bam.bai"), 
          emit: sanitized_bam
    path "versions.yml", emit: versions 

    script:
    def out_bam  = "${bam.baseName}.nochr.bam"
    
    """
    set -euo pipefail

    # Use multi-threading to remove 'chr' prefix from header
    if samtools view -H "${bam}" | grep "SN:chr" > /dev/null; then
        samtools view -@ ${task.cpus} -H "${bam}" \\
            | sed 's/SN:chr/SN:/g' \\
            | samtools reheader - "${bam}" > "${out_bam}"
    else
        # Copy original file if no changes needed
        cp "${bam}" "${out_bam}"
    fi

    # Generate BAM index using multiple threads
    samtools index -@ ${task.cpus} "${out_bam}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def out_bam  = "${bam.baseName}.nochr.bam"

    """
    touch ${out_bam}
    touch ${out_bam}.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | awk '{print \$2}')
    END_VERSIONS
    """
}