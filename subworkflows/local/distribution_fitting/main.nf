include { DISTFIT_HISTONE_TRAIN  } from '../../../modules/local/distfit_histone_train/main.nf'
include { DISTFIT_HISTONE_ASSESS } from '../../../modules/local/distfit_histone_assess/main.nf'

workflow DISTRIBUTION_FITTING {

    take:
    ch_input_raw  
    val_dists
    original_samplesheet

    main:
    ch_versions = Channel.empty()

    // Extract histone data
    ch_histone_data = ch_input_raw
        .map { it -> tuple(it[0], it[1]) }
        .filter { it[1] }

    ch_histone_inputs = ch_histone_data
        .flatMap { meta, h ->
            def header = h.withReader { it.readLine() }
            def marks_list = header.trim().split(/\s+/).toList().drop(3)
            marks_list.collectMany { mark ->
                val_dists.collect { dist -> tuple(meta, h, mark, dist) }
            }
        }

    // Call Training module
    DISTFIT_HISTONE_TRAIN(ch_histone_inputs)
    ch_versions = ch_versions.mix(DISTFIT_HISTONE_TRAIN.out.versions)

    // Prepare assess input
    ch_h_assess_input = ch_histone_data
        .map { meta, h -> tuple(meta.id, tuple(meta, h)) }
        .cross(DISTFIT_HISTONE_TRAIN.out.model.map { m, f -> tuple(m.id, f) }.groupTuple())
        .map { input, models -> tuple(input[1][0], input[1][1], models[1]) }
        
    // Call Assess process
    DISTFIT_HISTONE_ASSESS(ch_h_assess_input, original_samplesheet)
    ch_versions = ch_versions.mix(DISTFIT_HISTONE_ASSESS.out.versions)

    // Combine individual sample CSVs into one final master samplesheet
    DISTFIT_HISTONE_ASSESS.out.updated_samplesheet
        .map { meta, csv -> csv }
        .collectFile(
            name: 'optimal_samplesheet.csv', 
            storeDir: "${params.outdir}/pipeline_info", 
            keepHeader: true, 
            skip: 1
        )

    emit:
    updated_samplesheet = DISTFIT_HISTONE_ASSESS.out.updated_samplesheet 
    versions            = ch_versions
}