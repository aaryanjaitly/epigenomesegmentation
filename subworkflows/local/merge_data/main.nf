include { MERGE_COUNTS }        from '../../../modules/local/merge_counts/main.nf'
include { SPLIT_MERGED_COUNTS } from '../../../modules/local/split_merged_counts/main.nf'

workflow MERGE_DATA {
    
    take:
    ch_histone_data
    ch_methyl_data
    ch_bins
    ref_file

    main:
    ch_versions = Channel.empty()

    // Group methylation files by sample ID
    ch_methyl_data
        .map { meta, file ->
            def meta_grouped = meta.subMap('id')
            return [ meta_grouped, meta, file ]
        }
        .groupTuple()
        .set { ch_grouped_methyl }

    // Group histone files by sample ID
    ch_histone_data
        .map { meta, file ->
            def meta_grouped = meta.subMap('id')
            return [ meta_grouped, meta, file ]
        }
        .groupTuple()
        .set { ch_grouped_histone }

    // Join methylation and histone channels on sample ID
    ch_integrated_input = ch_grouped_methyl.join(ch_grouped_histone)

    // Map modality-specific data to genomic bins
    MERGE_COUNTS(
        ch_integrated_input,
        params.genome,
        ref_file,
        ch_bins
    )
    ch_versions = ch_versions.mix(MERGE_COUNTS.out.versions)

    // Separate merged data into histone and methylation model inputs
    SPLIT_MERGED_COUNTS(
        MERGE_COUNTS.out.merged_counts
    )
    ch_versions = ch_versions.mix(SPLIT_MERGED_COUNTS.out.versions)

    emit:
    split_counts = SPLIT_MERGED_COUNTS.out.inputs
    versions     = ch_versions
}