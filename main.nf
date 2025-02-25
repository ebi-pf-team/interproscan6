nextflow.enable.dsl=2

include { INIT_PIPELINE                 } from "./interproscan/subworkflows/init"
include { SCAN_SEQUENCES                } from "./interproscan/subworkflows/scan"

include { POPULATE_SEQ_DATABASE;
          UPDATE_ORFS;
          BUILD_BATCHES;
          INDEX_BATCHES                 } from "./interproscan/modules/prepare_sequences"
include { ESL_TRANSLATE                 } from "./interproscan/modules/esl_translate"
include { LOOKUP_MATCHES                } from "./interproscan/modules/lookup"
include { XREFS                         } from "./interproscan/modules/xrefs"
include { REPRESENTATIVE_DOMAINS        } from "./interproscan/modules/representative_domains"
include { WRITE_JSON_OUTPUT             } from "./interproscan/modules/output/json"
include { WRITE_TSV_OUTPUT              } from "./interproscan/modules/output/tsv"
include { WRITE_XML_OUTPUT              } from "./interproscan/modules/output/xml"


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
    formats         = INIT_PIPELINE.out.formats.val
    apps            = INIT_PIPELINE.out.apps.val
    signalpMode     = INIT_PIPELINE.out.signalpMode.val

    if (params.nucleic) {
        // Store the input seqs in the internal ips6 seq db
        db_path = POPULATE_SEQ_DATABASE(fasta_file, params.nucleic)

        // Chunk input file in smaller files for translation
        fasta_file
            .splitFasta( by: params.batchSize, file: true )
            .set { ch_fasta }

        // Translate DNA/RNA sequences to protein sequences
        ESL_TRANSLATE(ch_fasta)

        // Store sequences in the sequence database
        UPDATE_ORFS(ESL_TRANSLATE.out, db_path)

        // Wait for all UPDATE_ORFS to complete before proceeding
        UPDATE_ORFS.out.collect().map { db_path }.set { ready_db_path }
    } else {
        // Store the input seqs in the internal ips6 seq db
        ready_db_path = POPULATE_SEQ_DATABASE(fasta_file, params.nucleic)
    }
    // Build batches of unique protein seqs for the analysis
    BUILD_BATCHES(ready_db_path, params.batchSize)

    // [fasta, fasta, fasta] --> [[index, fasta], [index, fasta], [index, fasta]] - to help gather matches for each batch
    ch_seqs = INDEX_BATCHES(BUILD_BATCHES.out).flatMap { it } // flatMap so tuples are emitted one at a time

    matchResults = Channel.empty()
    if (params.disablePrecalc) {
        SCAN_SEQUENCES(
            ch_seqs,
            apps,
            params.appsConfig,
            data_dir
        )
        matchResults = SCAN_SEQUENCES.out
    } else {
        LOOKUP_MATCHES(
            ch_seqs,
            apps,
            params.lookupService.apiChunkSize,
            params.lookupService.lookupHost,
            params.lookupService.maxRetries
        )

        SCAN_SEQUENCES(
            LOOKUP_MATCHES.out[1],
            apps,
            params.appsConfig,
            data_dir
        )

        def expandedScan = SCAN_SEQUENCES.out.flatMap { scan ->
            scan[1].collect { path -> [scan[0], path] }
        }

        def combined = LOOKUP_MATCHES.out[0].concat(expandedScan)
        matchResults = combined.groupTuple()
    }
    // matchResults = [[meta, [member.json, member.json, member.json]]

    /* XREFS:
    Aggregate matches across all members for each sequence --> single JSON with all matches for the batch
    Add signature and entry desc and names
    Add PAINT annotations (if panther is enabled)
    Add go terms (if enabled)
    Add pathways (if enabled)
    */
    XREFS(
        matchResults,
        apps,
        data_dir,
        params.xRefsConfig.entries,
        params.xRefsConfig.goterms,
        params.xRefsConfig.pathways,
        params.goterms,
        params.pathways,
        "${data_dir}/${params.appsConfig.paint}"
    )

    REPRESENTATIVE_DOMAINS(XREFS.out)
    // Collect all JSON files into a single channel so we don't have cocurrent writing to the output files
    ch_results = REPRESENTATIVE_DOMAINS.out
        .map { meta, json -> json }
        .collect()

    def fileName = params.input.split('/').last()
    def outFileName = "${params.outdir}/${fileName}"

// //     if (formats.contains("JSON")) {
// //         WRITE_JSON_OUTPUT(Rch_results, "${outFileName}", params.nucleic, workflow.manifest.version)
// //     }
    if (formats.contains("TSV")) {
        WRITE_TSV_OUTPUT(ch_results, "${outFileName}", db_path, params.nucleic)
    }
// //     if (formats.contains("XML")) {
// //         WRITE_XML_OUTPUT(ch_results, "${outFileName}", params.nucleic, workflow.manifest.version)
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