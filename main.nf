nextflow.enable.dsl=2

include { INIT_PIPELINE                 } from "./interproscan/subworkflows/init"
include { SCAN_SEQUENCES                } from "./interproscan/subworkflows/scan"

include { POPULATE_SEQ_DATABASE;
          UPDATE_ORFS;
          BUILD_BATCHES;
          INDEX_FASTA_FILES             } from "./interproscan/modules/prepare_sequences"
include { ESL_TRANSLATE                 } from "./interproscan/modules/esl_translate"
include { LOOKUP_MATCHES                } from "./interproscan/modules/lookup"
include { XREFS                         } from "./interproscan/modules/xrefs"
include { AGGREGATE_MATCHES             } from "./interproscan/modules/aggregate_matches"
include { REPRESENTATIVE_DOMAINS        } from "./interproscan/modules/representative_domains"
include { WRITE_JSON_OUTPUT             } from "./interproscan/modules/output/json"
include { WRITE_TSV_OUTPUT              } from "./interproscan/modules/output/tsv"
include { WRITE_XML_OUTPUT } from "./interproscan/modules/output/xml"


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

    if (params.nucleic) {
        // Store the input seqs in the internal ips6 seq db
        POPULATE_SEQ_DATABASE(fasta_file, params.nucleic)

        // Chunk input file in smaller files for translation
        fasta_file
            .splitFasta( by: params.batchSize, file: true )
            .set { ch_fasta }

        // Translate DNA/RNA sequences to protein sequences
        ch_translated = ESL_TRANSLATE(ch_fasta)

        // Store sequences in the sequence database
        UPDATE_ORFS(ch_translated, POPULATE_SEQ_DATABASE.out)

        // Build batches of unique protein seqs for the analysis
        fasta_files = BUILD_BATCHES(UPDATE_ORFS.out, params.batchSize)
    } else {
        // Store the input seqs in the internal ips6 seq db
        POPULATE_SEQ_DATABASE(fasta_file, params.nucleic)

        // Build batches of unique protein seqs for the analysis
        fasta_files = BUILD_BATCHES(POPULATE_SEQ_DATABASE.out, params.batchSize)
    }

    // Index fasta files [fasta, fasta, fasta] --> [[index, fasta], [index, fasta]]
    // This aids aggregating matches from all member databases for each batch file
    ch_seqs = INDEX_FASTA_FILES(BUILD_BATCHES.out)
    ch_seqs.view()
//     ch_seqs.view()
//     LOOKUP_MATCHES(
//         ch_seqs,
//         apps,
//         params.lookupService.apiChunkSize,
//         params.lookupService.lookupHost,
//         params.lookupService.maxRetries
//     )

//
//     matchResults = Channel.empty()
//     if (params.disablePrecalc) {
//         SCAN_SEQUENCES(
//             BUILD_BATCHES.out,
//             apps,
//             params.appsConfig,
//             data_dir
//         )
//         matchResults = SCAN_SEQUENCES.out
//     } else {
//         LOOKUP_MATCHES(
//             BUILD_BATCHES.out,
//             apps,
//             params.lookupService.apiChunkSize,
//             params.lookupService.lookupHost,
//             params.lookupService.maxRetries
//         )
//
//         SCAN_SEQUENCES(
//             LOOKUP_MATCHES.out[1],  // [index, fasta of seqs not in the MLS]
//             apps,
//             params.appsConfig,
//             data_dir
//         )
//
//         def expandedScan = SCAN_SEQUENCES.out.flatMap { scan ->
//             scan[1].collect { path -> [scan[0], path] }
//         }
//
//         def combined = LOOKUP_MATCHES.out[0].concat(expandedScan)
//         matchResults = combined.groupTuple()
//     }
//
//     /* XREFS:
//     Add signature and entry desc and names
//     Add PAINT annotations (if panther is enabled)
//     Add go terms (if enabled)
//     Add pathways (if enabled)
//     */
//     XREFS(
//         matchResults,
//         apps,
//         data_dir,
//         params.xRefsConfig.entries,
//         params.xRefsConfig.goterms,
//         params.xRefsConfig.pathways,
//         params.goterms,
//         params.pathways,
//         "${data_dir}/${params.appsConfig.paint}"
//     )
//
//     ch_seqs.join(XREFS.out, by: 0)
//     .map { batchnumber, fasta, sequences, matches ->
//         [batchnumber, sequences, matches]
//     }.set { ch_seq_matches }
//
//     AGGREGATE_MATCHES(ch_seq_matches, params.nucleic)
//
//     REPRESENTATIVE_DOMAINS(AGGREGATE_MATCHES.out)
//
// //     Channel.from(params.formats.toLowerCase().split(','))
// //     .set { ch_format }
// //
// //     def formats = params.formats.toUpperCase().split(',') as Set
// //     def fileName = params.input.split('/').last()
// //     def outFileName = "${params.outdir}/${fileName}"
// //     if (formats.contains("JSON")) {
// //         WRITE_JSON_OUTPUT(REPRESENTATIVE_DOMAINS.out, "${outFileName}", params.nucleic, workflow.manifest.version)
// //     }
// //     if (formats.contains("TSV")) {
// //         WRITE_TSV_OUTPUT(REPRESENTATIVE_DOMAINS.out, "${outFileName}", params.nucleic)
// //     }
// //     if (formats.contains("XML")) {
// //         WRITE_XML_OUTPUT(REPRESENTATIVE_DOMAINS.out, "${outFileName}", params.nucleic, workflow.manifest.version)
// //     }
}

workflow.onComplete = {
    def input_file = file(params.input)
    def outputFileName = input_file.getName()
    def outputDir = params.outdir.endsWith('/') ? params.outdir[0..-2] : params.outdir

    println "InterProScan workflow completed successfully: $workflow.success."
    println "Any results are located at ${outputDir}/${outputFileName}.ips6.*"
    println "Duration: $workflow.duration"
}