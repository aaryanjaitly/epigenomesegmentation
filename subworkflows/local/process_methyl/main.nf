include { GENERATE_COUNT_MATRIX_BED } from '../../../modules/local/generate_count_matrix_bed/main.nf'

workflow PROCESS_METHYL {
    
    take:
    ch_input_bed

    main:
    ch_versions = Channel.empty()

    // Generate methylation count matrices from input BED files
    GENERATE_COUNT_MATRIX_BED(
        ch_input_bed
    )
    ch_versions = ch_versions.mix(GENERATE_COUNT_MATRIX_BED.out.versions)

    emit:
    counts   = GENERATE_COUNT_MATRIX_BED.out.counts_bed
    versions = ch_versions
}