include { GENERATE_METHYLATION_BINS }   from '../../../modules/local/generate_methylation_bins/main.nf'

workflow GENERATE_BINS {
    
    take:
    ch_methyl_data
    ch_bins

    main:
    ch_versions = Channel.empty()

    ch_methyl_grouped = ch_methyl_data
        .map { meta, tab ->
            // Pull 'modality' straight out of the meta map
            def mod = meta.modality 
            
            // Return the new 3-item tuple
            [ meta, mod, tab ]
        }
        .groupTuple(by: 0)

        
    GENERATE_METHYLATION_BINS(
        ch_methyl_grouped.map { meta, mods, tabs ->
            [ 
              meta, 
              mods, 
              tabs 
            ]
        },
        ch_bins
    )

    ch_versions = ch_versions.mix(GENERATE_METHYLATION_BINS.out.versions)

    emit:
    split_counts = GENERATE_METHYLATION_BINS.out.merged_counts
    versions     = ch_versions
}