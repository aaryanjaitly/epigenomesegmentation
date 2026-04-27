include { DNA_PREPARE } from '../../../modules/local/dna_prepare/main.nf'
include { DNA_TRAIN   } from '../../../modules/local/dna_train/main.nf'
include { DNA_DECODE  } from '../../../modules/local/dna_decode/main.nf'
include { DNA_REPORT  } from '../../../modules/local/dna_report/main.nf'

workflow MODEL_TRAINING_DNA {

    take:
    ch_input

    main:
    ch_versions = Channel.empty()

    ch_inputs_split = ch_input.flatMap { meta, meth, states ->

        def state_list = (states instanceof String)
            ? states.split(',')
            : states

        state_list.collect { state ->
            def new_meta = meta.clone()
            new_meta.id    = "${meta.id}_s${state}"
            new_meta.state = state

            [ new_meta, meth, state ]
        }
    }

    DNA_PREPARE(ch_inputs_split)
    ch_versions = ch_versions.mix(DNA_PREPARE.out.versions)

    ch_train_inputs = DNA_PREPARE.out.config
        .join(DNA_PREPARE.out.train_counts)
        .join(DNA_PREPARE.out.train_regions)

    DNA_TRAIN(ch_train_inputs)
    ch_versions = ch_versions.mix(DNA_TRAIN.out.versions)

    ch_decode_inputs = DNA_PREPARE.out.config
        .join(DNA_TRAIN.out.model)
        .join(DNA_PREPARE.out.counts)
        .join(DNA_PREPARE.out.train_regions)
        
    DNA_DECODE(ch_decode_inputs)
    ch_versions = ch_versions.mix(DNA_DECODE.out.versions)

    ch_report_inputs = DNA_PREPARE.out.config
        .join(DNA_TRAIN.out.model)
        .join(DNA_DECODE.out.bed)
        .join(DNA_DECODE.out.seg_txt)

    DNA_REPORT(ch_report_inputs)
    ch_versions = ch_versions.mix(DNA_REPORT.out.versions)

    emit:
    versions = ch_versions
}