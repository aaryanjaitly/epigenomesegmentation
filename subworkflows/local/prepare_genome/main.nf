include { GET_GENOME_CHROM_SIZE } from '../../../modules/local/get_genome_chrom_size/main.nf'
include { GENERATE_GENOME_BINS } from '../../../modules/local/generate_genome_bins/main.nf'

workflow PREPARE_GENOME {
    main:
    ch_versions = Channel.empty()

    // Download and filter chromosome sizes for the specified genome
    GET_GENOME_CHROM_SIZE(
        params.genome
    )
    ch_versions = ch_versions.mix(GET_GENOME_CHROM_SIZE.out.versions)

    // Create fixed-size genomic windows from chromosome sizes
    GENERATE_GENOME_BINS(
        GET_GENOME_CHROM_SIZE.out.genome_chrom_size_file,
        params.genome,
        params.binsize
    )
    ch_versions = ch_versions.mix(GENERATE_GENOME_BINS.out.versions)

    emit:
    chrom_sizes = GET_GENOME_CHROM_SIZE.out.genome_chrom_size_file
    bins        = GENERATE_GENOME_BINS.out.bins_file
    versions    = ch_versions
}