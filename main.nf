nextflow.enable.dsl=2

include { INIT_PIPELINE                 } from "./interproscan/subworkflows/init"
include { SCAN_SEQUENCES                } from "./interproscan/subworkflows/scan"
include { ESL_TRANSLATE                 } from "./interproscan/modules/esl_translate"
include { PREPARE_NUCLEIC_SEQUENCES     } from "./interproscan/modules/prepare_sequences"
include { PREPARE_PROTEIN_SEQUENCES     } from "./interproscan/modules/prepare_sequences"

workflow {
    println "# ${workflow.manifest.name} ${workflow.manifest.version}"
    println "# ${workflow.manifest.description}\n"

    if (params.keySet().any { it.equalsIgnoreCase("help") }) {
        InterProScan.printHelp(params.appsConfig)
        exit 0
    }

    INIT_PIPELINE()

    fasta_file      = Channel.fromPath(INIT_PIPELINE.out.fasta.val)
    data_dir        = INIT_PIPELINE.out.datadir.val
    outut_dir       = INIT_PIPELINE.out.outdir.val
    apps            = INIT_PIPELINE.out.apps.val
    signalpMode     = INIT_PIPELINE.out.signalpMode.val

    // Chunk input file in smaller files
    fasta_file
        .splitFasta( by: params.batchSize, file: true )
        .set { ch_fasta }

    if (params.nucleic) {
        // Translate DNA/RNA sequences to protein sequences
        ch_translated = ESL_TRANSLATE(ch_fasta)

        // Split again
        ch_translated
            .splitFasta( by: params.batchSize, file: true )
            .map { split_pt_file, orig_nt_file -> tuple ( orig_nt_file, split_pt_file ) }
            .set { ch_translated_split }

        // Store sequences as JSON objects
        ch_seqs = PREPARE_NUCLEIC_SEQUENCES(ch_translated_split)
    } else {
        // Store sequences as JSON objects
        ch_seqs = PREPARE_PROTEIN_SEQUENCES(ch_fasta)
    }

    // ch_seqs
    //     .map { index, fasta, json -> tuple( index, fasta ) }
    //     .set { ch_fasta }

    // ch_seqs
    //     .map { index, fasta, json -> json }
    //     .set { ch_json }

    // TODO: add new match lookup
    
    SCAN_SEQUENCES(
        ch_seqs,
        apps,
        params.appsConfig,
        data_dir,
        signalpMode
    )

    SCAN_SEQUENCES.out.view()


    // disable_precalc = params.disable_precalc
    // sequences_to_analyse = null
    // parsed_matches = Channel.empty()
    // if (!disable_precalc) {
    //     log.info "Using precalculated match lookup service"
    //     SEQUENCE_PRECALC(PARSE_SEQUENCE.out, applications, false)  // final: bool to indicate not a unit test
    //     sequences_to_analyse = SEQUENCE_PRECALC.out.sequences_to_analyse
    //     parsed_matches = SEQUENCE_PRECALC.out.parsed_matches
    // }

    // if (parsed_matches.collect() == null) {
    //         // cases in which the lookup check ran successfully but lookup matches not
    //         disable_precalc = true
    //         log.info "ERROR: unable to connect to match lookup service. Max retries reached. Running analysis locally..."
    // }

    // analysis_result = Channel.empty()
    // if (disable_precalc || sequences_to_analyse) {
    //     log.info "Running sequence analysis"
    //     if (sequences_to_analyse && !disable_precalc) {
    //         fasta_to_runner = sequences_to_analyse
    //     }
    //     else {
    //         if (params.nucleic) {
    //             fasta_to_runner = orfs_fasta
    //         }
    //         else {
    //             fasta_to_runner = ch_fasta
    //         }
    //     }
    //     parsed_analysis = SEQUENCE_ANALYSIS(
    //         fasta_to_runner,
    //         applications,
    //         dataDirPath,
    //         params.signalp_mode
    //     )
    // }

    // AGGREGATE_PARSED_SEQS(PARSE_SEQUENCE.out.collect())

    // all_results = parsed_matches.concat(parsed_analysis)
    // AGGREGATE_RESULTS(all_results)

    // /* XREFS:
    // Add signature and entry desc and names
    // Add PAINT annotations (if panther is enabled)
    // Add go terms (if enabled)
    // Add pathways (if enabled)
    // */
    // XREFS(AGGREGATE_RESULTS.out, applications, dataDirPath)

    // REPRESENTATIVE_DOMAINS(XREFS.out.collect())

    // Channel.from(params.formats.toLowerCase().split(','))
    // .set { ch_format }

    // WRITE_RESULTS(
    //     input_file.getName(),
    //     AGGREGATE_PARSED_SEQS.out,
    //     REPRESENTATIVE_DOMAINS.out.collect(),
    //     ch_format,
    //     params.outdir,
    //     params.ipscn_version,
    //     params.nucleic
    // )
}

// workflow.onComplete = {
//     def input_file = file(params.input)
//     def outputFileName = input_file.getName()
//     def outputDir = params.outdir.endsWith('/') ? params.outdir[0..-2] : params.outdir

//     println "InterProScan workflow completed successfully: $workflow.success."
//     println "Any results are located at ${outputDir}/${outputFileName}.ips6.*"
//     println "Duration: $workflow.duration"
// }
