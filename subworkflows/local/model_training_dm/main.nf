include { DM_PREPARE} from '../../../modules/local/dm_prepare/main.nf'
include { DM_TRAIN   } from '../../../modules/local/dm_train/main.nf'
include { DM_DECODE  } from '../../../modules/local/dm_decode/main.nf'
include { DM_REPORT  } from '../../../modules/local/dm_report/main.nf'

workflow MODEL_TRAINING_DM {

    take:
    ch_input       
     
    main:
    ch_versions = Channel.empty()

    // Create unique channel entries for each state number
    ch_dm_prepare_in = ch_input
        .flatMap { meta, histone, meth ->
            def state_list = (params.states instanceof String) ? params.states.split(',') : params.states
            
            state_list.collect { state ->
                def new_meta = meta.clone()
                new_meta.id = "${meta.id}_s${state}"
                
                [ new_meta, histone, meth, state ]
            }
        }

    // Generate configuration and count matrices
    DM_PREPARE(
        ch_dm_prepare_in,    
        params.chr_parameter_estimation  
    )
    ch_versions = ch_versions.mix(DM_PREPARE.out.versions)

    // Collect training data components
    ch_train_inputs = DM_PREPARE.out.config
        .join(DM_PREPARE.out.train_counts)
        .join(DM_PREPARE.out.regions)
        .join(DM_PREPARE.out.train_meth)

    // Execute HMM training
    DM_TRAIN(ch_train_inputs)
    ch_versions = ch_versions.mix(DM_TRAIN.out.versions)

    // Combine model and config for Viterbi decoding
    ch_decode_inputs = DM_PREPARE.out.config
        .join(DM_TRAIN.out.model)

    // Perform chromatin segmentation
    DM_DECODE(ch_decode_inputs)
    ch_versions = ch_versions.mix(DM_DECODE.out.versions)

    // Assemble all results for final visualization and reporting
    ch_report_inputs = DM_PREPARE.out.config
        .join(DM_TRAIN.out.model)
        .join(DM_DECODE.out.bed)
        .join(DM_DECODE.out.seg_txt)

    // Generate plots and HTML summary
    DM_REPORT(ch_report_inputs)
    ch_versions = ch_versions.mix(DM_REPORT.out.versions)

    emit:
    versions = ch_versions
}