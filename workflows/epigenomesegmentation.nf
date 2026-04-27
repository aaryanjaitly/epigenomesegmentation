/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_epigenomesegmentation_pipeline'

/* LOCAL SUBWORKFLOWS */
include { PREPARE_GENOME       } from '../subworkflows/local/prepare_genome/main.nf'
include { PROCESS_METHYL       } from '../subworkflows/local/process_methyl/main.nf'
include { PROCESS_HISTONES     } from '../subworkflows/local/process_histones/main.nf'
include { MERGE_DATA           } from '../subworkflows/local/merge_data/main.nf'
include { GENERATE_BINS        } from '../subworkflows/local/generate_bins/main.nf'
include { MODEL_TRAINING_STD   } from '../subworkflows/local/model_training_std/main.nf'
include { MODEL_TRAINING_DM    } from '../subworkflows/local/model_training_dm/main.nf'
include { MODEL_TRAINING_DNA   } from '../subworkflows/local/model_training_dna/main.nf'
include { DISTRIBUTION_FITTING } from '../subworkflows/local/distribution_fitting/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EPIGENOMESEGMENTATION {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    
    main:

    ch_versions = channel.empty()

    // =========================================================
    // START OF INJECTED EPISEGMIX LOGIC
    // =========================================================
    ch_counts        = Channel.empty()

    // ---------------------------------------------------------
    // DETERMINE PIPELINE BEHAVIOR FLAGS
    // ---------------------------------------------------------
    def is_dna      = params.DNA
    def is_counts   = params.counts

    // ---------------------------------------------------------
    // FORK SAMPLESHEET TO PREVENT CHANNEL CONSUMPTION HANGS
    // ---------------------------------------------------------
    ch_samplesheet.multiMap { it ->
        for_dists: it
        for_branching: it
    }.set { ch_forked_samplesheet }

    // ---------------------------------------------------------
    // INITIAL DISTRIBUTION MAP
    // ---------------------------------------------------------
    ch_distributions = ch_forked_samplesheet.for_dists
        .map { tuple -> 
            def meta = tuple[0]
            def dist = meta.distribution ?: meta.distributions
            [ meta.id, meta.epigenetic_mark, dist ] 
        }
        .groupTuple(by: 0)
        .map { tuple ->
            def id    = tuple[0]
            def marks = tuple[1]
            def dists = tuple[2]
            def dist_map = [marks, dists].transpose().collectEntries { k, v -> [k, v] }
            return [ id, dist_map ]
        }

    // ---------------------------------------------------------
    // STEP 1: PREPARE GENOME
    // ---------------------------------------------------------
    PREPARE_GENOME()
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    // ---------------------------------------------------------
    // STEP 2: SPLIT INPUT
    // ---------------------------------------------------------
    ch_forked_samplesheet.for_branching
        .branch { it ->
            bam: it[1].toString().toLowerCase().endsWith('.bam') || it[1].toString().toLowerCase().endsWith('.bam.gz')
            bed: it[1].toString().toLowerCase().endsWith('.bed') || it[1].toString().toLowerCase().endsWith('.bed.gz')
            other: true
        }.set { ch_input_branched }

    // ---------------------------------------------------------
    // STEP 3: HISTONES (Skip if pure DNA mode)
    // ---------------------------------------------------------
    if (!is_dna) {
        PROCESS_HISTONES(ch_input_branched.bam, PREPARE_GENOME.out.chrom_sizes)
        ch_versions = ch_versions.mix(PROCESS_HISTONES.out.versions)
        ch_counts   = ch_counts.mix(PROCESS_HISTONES.out.counts)
    }
    
    // ---------------------------------------------------------
    // STEP 4: METHYLATION
    // ---------------------------------------------------------
    if (is_dna || params.merge) {
        PROCESS_METHYL(ch_input_branched.bed)
        ch_versions = ch_versions.mix(PROCESS_METHYL.out.versions)
        ch_counts   = ch_counts.mix(PROCESS_METHYL.out.counts)
    }

    // ---------------------------------------------------------
    // CONDITIONAL EXECUTION (SKIPPED IN 'COUNTS' MODE)
    // ---------------------------------------------------------
    if (!is_counts) {

        // ---------------------------------------------------------
        // STEP 5: MERGING / BINNING
        // ---------------------------------------------------------
        ch_model_input = Channel.empty()

        if (is_dna) {
            GENERATE_BINS(PROCESS_METHYL.out.counts, PREPARE_GENOME.out.bins)
            ch_model_input = GENERATE_BINS.out.split_counts
            ch_versions    = ch_versions.mix(GENERATE_BINS.out.versions)

        } else if (params.merge) {
            MERGE_DATA(
                PROCESS_HISTONES.out.counts,
                PROCESS_METHYL.out.counts,
                PREPARE_GENOME.out.bins,
                PREPARE_GENOME.out.chrom_sizes
            )
            ch_model_input = MERGE_DATA.out.split_counts
            ch_versions    = ch_versions.mix(MERGE_DATA.out.versions)

        } else {
            ch_model_input = PROCESS_HISTONES.out.counts
                .map { tuple -> [ tuple[0], tuple[1], [] ] }
        }

        // ---------------------------------------------------------
        // STEP 6: INITIAL METADATA INJECTION
        // ---------------------------------------------------------
        ch_ready_for_analysis = ch_model_input
            .map { tuple -> [ tuple[0].id, tuple ] } 
            .join(ch_distributions)                  
            .map { tuple ->
                def id             = tuple[0]
                def original_tuple = tuple[1]
                def dist_map       = tuple[2]
                def new_meta       = original_tuple[0] + [ distributions: dist_map ]    
                return [ new_meta ] + original_tuple.drop(1) 
            }

        // ---------------------------------------------------------
        // STEP 7: SEQUENTIAL MODEL ROUTING
        // ---------------------------------------------------------
        // Default to the base analysis input
        def ch_final_training_input = ch_ready_for_analysis

        // PHASE A: FITTING
        if (params.fitting) {
            def dist_list = params.distributions ? params.distributions.split(',').collect{ it.trim() } : []
            
            DISTRIBUTION_FITTING(
                ch_ready_for_analysis.map { tuple -> [tuple[0], tuple[1]] },
                dist_list,
                Channel.fromPath(params.input).first() 
            )
            ch_versions = ch_versions.mix(DISTRIBUTION_FITTING.out.versions)

            // CHAINING LOGIC: If segmentation is requested, parse the new CSV and update the input
            if (params.best_fit_segmentation) {
                
                // Parse the new distributions from the CSV output
                ch_updated_dist = DISTRIBUTION_FITTING.out.updated_samplesheet
                    .splitCsv(header: true)
                    .map { tuple -> 
                        def meta = tuple[0]
                        def row  = tuple[1]
                        [ meta.id, row.epigenetic_mark, row.distribution ] 
                    }
                    .groupTuple(by: 0)
                    .map { tuple ->
                        def id    = tuple[0]
                        def marks = tuple[1]
                        def dists = tuple[2]
                        def dist_map = [marks, dists].transpose().collectEntries { k, v -> [k, v] }
                        return [ id, dist_map ]
                    }

                // Inject the updated distributions back into the main data channel
                ch_final_training_input = ch_ready_for_analysis
                    .map { tuple -> [ tuple[0].id, tuple ] } 
                    .join(ch_updated_dist)                  
                    .map { tuple ->
                        def id             = tuple[0]
                        def original_tuple = tuple[1]
                        def dist_map       = tuple[2]
                        def new_meta       = original_tuple[0] + [ distributions: dist_map ]    
                        return [ new_meta ] + original_tuple.drop(1) 
                    }
            }
        }

        // PHASE B: SEGMENTATION / DURATION / STANDARD TRAINING
        // Runs if it's a standard run, OR if it's a fitting run that explicitly requested downstream segmentation
        if (!params.fitting || params.best_fit_segmentation) {
            
            if (params.duration) {
                MODEL_TRAINING_DM(ch_final_training_input)
                ch_versions = ch_versions.mix(MODEL_TRAINING_DM.out.versions)

            } else if (is_dna) {
                ch_dna_model_input = ch_final_training_input.map { tuple -> 
                    [ tuple[0], tuple[1], tuple.size() > 2 ? tuple[2] : params.states ] 
                }
                MODEL_TRAINING_DNA(ch_dna_model_input)
                ch_versions = ch_versions.mix(MODEL_TRAINING_DNA.out.versions)

            } else {
                MODEL_TRAINING_STD(ch_final_training_input)
                ch_versions = ch_versions.mix(MODEL_TRAINING_STD.out.versions)
            }
        }
    }
    // =========================================================
    // END OF INJECTED EPISEGMIX LOGIC
    // =========================================================

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { tuple ->
            def process = tuple[0]
            def tool    = tuple[1]
            def version = tuple[2]
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { tuple ->
            def process       = tuple[0]
            def tool_versions = tuple[1]
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'epigenomesegmentation_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
