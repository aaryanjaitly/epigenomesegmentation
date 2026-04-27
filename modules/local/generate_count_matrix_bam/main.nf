process GENERATE_COUNT_MATRIX_BAM {
    tag "Counts: ${meta.id}"
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://aaryanjaitly/episegmix_counts_container:latest' :
        'aaryanjaitly/episegmix_counts_container:latest' }"

    input:
    tuple val(meta), 
          val(meta_list), 
          path(bams), 
          path(bais)
    path reference_file 
    val  genome
    
    output:
    tuple val(meta), 
          path("${meta.id}_Histones/*_refined_counts.txt"), 
          emit: counts
    path "versions.yml", emit: versions

    script:
    def sample_id = meta.id
    def out_dir   = "${sample_id}_Histones"
    def pe_lines  = []
    def se_lines  = []

    meta_list.eachWithIndex { m, index ->
        def mark     = m.epigenetic_mark
        def is_pe    = m.paired_end
        def bam_file = bams instanceof List ? bams[index].name : bams.name
        def line     = "${mark}\t${bam_file}"
        
        if (is_pe) { 
            pe_lines.add(line) 
        } else { 
            se_lines.add(line) 
        }
    }

    def pe_content = pe_lines.join('\n')
    def se_content = se_lines.join('\n')

    """
    set -euo pipefail
    mkdir -p "${out_dir}"
    
    PE_OUT="${out_dir}/${sample_id}_PE_${genome}_refined_counts.txt"
    SE_OUT="${out_dir}/${sample_id}_SE_${genome}_refined_counts.txt"
    FINAL_OUT="${out_dir}/${sample_id}_${genome}_refined_counts.txt"

    # Create a base reference file that won't be touched by the scripts
    awk -v OFS="\t" '{print \$1, 1, \$2}' "${reference_file}" > new_reference_base

    if [[ -n "${pe_content}" ]]; then
        echo "${pe_content}" > "PE_info.txt"
        
        # Create a specific copy for the PE run
        cp new_reference_base new_reference_pe
        
        bash counts.sh \\
            -t "PE_info.txt" \\
            -o "${out_dir}" \\
            -c ${task.cpus} \\
            -p midpoint \\
            -f "${sample_id}_PE" \\
            -g "${genome}" \\
            -r ./new_reference_pe \\
            -b ${params.binsize}
    fi

    if [[ -n "${se_content}" ]]; then
        echo "${se_content}" > "SE_info.txt"
        
        # Create a specific copy for the SE run
        cp new_reference_base new_reference_se
        
        bash counts.sh \\
            -t "SE_info.txt" \\
            -o "${out_dir}" \\
            -c ${task.cpus} \\
            -p ignore \\
            -f "${sample_id}_SE" \\
            -g "${genome}" \\
            -r ./new_reference_se \\
            -b ${params.binsize}
    fi

    if [[ -f "\$PE_OUT" && -f "\$SE_OUT" ]]; then
        paste "\$PE_OUT" <(cut -f4- "\$SE_OUT") > "\$FINAL_OUT"
        rm "\$PE_OUT" "\$SE_OUT"
    elif [[ -f "\$PE_OUT" ]]; then
        mv "\$PE_OUT" "\$FINAL_OUT"
    elif [[ -f "\$SE_OUT" ]]; then
        mv "\$SE_OUT" "\$FINAL_OUT"
    else
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n 1 | sed 's/.*version //g' | awk '{print \$1}')
        awk: \$(awk --version | head -n 1 | awk '{print \$3}')
END_VERSIONS
    """

    stub:
    def sample_id = meta.id
    def out_dir   = "${sample_id}_Histones"

    """
    mkdir -p "${out_dir}"
    touch "${out_dir}/${sample_id}_${genome}_refined_counts.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n 1 | sed 's/.*version //g' | awk '{print \$1}')
        awk: \$(awk --version | head -n 1 | awk '{print \$3}')
END_VERSIONS
    """
}