include { CLEAN_AND_INDEX_BAM }       from '../../../modules/local/clean_and_index_bam/main.nf'
include { GENERATE_COUNT_MATRIX_BAM } from '../../../modules/local/generate_count_matrix_bam/main.nf'

workflow PROCESS_HISTONES {
    
    take:
    ch_input_bam
    chrom_sizes

    main:
    ch_versions = Channel.empty()
    
    // Clean headers and index input BAM files
    CLEAN_AND_INDEX_BAM(
        ch_input_bam
    )
    ch_versions = ch_versions.mix(CLEAN_AND_INDEX_BAM.out.versions)

    // Group processed BAMs into sample-specific tuples for matrix generation
    CLEAN_AND_INDEX_BAM.out.sanitized_bam
        .map { tuple ->
            def meta = tuple[0]
            def bam  = tuple[1]
            def bai  = tuple[2]
            def group_key = meta.subMap('id')
            return [ group_key, [ meta, bam, bai ] ]
        }
        .groupTuple()
        .map { tuple ->
            def group_key = tuple[0]
            def inputs    = tuple[1]
            def list_of_metas = inputs.collect { it[0] }
            def list_of_bams  = inputs.collect { it[1] }
            def list_of_bais  = inputs.collect { it[2] }

            return [ group_key, list_of_metas, list_of_bams, list_of_bais ]
        }
        .set { ch_grouped_bams }

    // Generate refined histone count matrices from grouped BAMs
    GENERATE_COUNT_MATRIX_BAM(
        ch_grouped_bams,
        chrom_sizes,
        params.genome
    )
    ch_versions = ch_versions.mix(GENERATE_COUNT_MATRIX_BAM.out.versions)
    
    emit:
    counts   = GENERATE_COUNT_MATRIX_BAM.out.counts
    versions = ch_versions
}