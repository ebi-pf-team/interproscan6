nextflow.enable.dsl=2
import groovy.json.JsonSlurper  // until selective downloads

include { INIT_PIPELINE                 } from "./interproscan/subworkflows/init"
include { PREPARE_SEQUENCES             } from "./interproscan/subworkflows/prepare_sequences"
include { SCAN_SEQUENCES                } from "./interproscan/subworkflows/scan"
include { INTERPRO                      } from "./interproscan/subworkflows/interpro"
include { OUTPUT                        } from "./interproscan/subworkflows/output"

include { LOOKUP_MATCHES                } from "./interproscan/modules/lookup" // will be migrated to a subworkflow when selective downloads is merged

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
    matchesApiUrl   = INIT_PIPELINE.out.matchesApiUrl.val

    PREPARE_SEQUENCES(fasta_file, apps)
    ch_seqs         = PREPARE_SEQUENCES.out.ch_seqs
    seq_db_path     = PREPARE_SEQUENCES.out.seq_db_path

    match_results = Channel.empty()
    if (matchesApiUrl != null) {
        LOOKUP_MATCHES(
            ch_seqs,
            apps,
            matchesApiUrl,
            params.lookupService.chunkSize,
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
        match_results = combined.groupTuple()
    } else {
        SCAN_SEQUENCES(
            ch_seqs,
            apps,
            params.appsConfig,
            data_dir
        )
        match_results = SCAN_SEQUENCES.out
    }
    // match_results format: [[meta, [member1.json, member2.json, ..., memberN.json]]

    /* INTERPRO:
    Aggregate matches across all members for each sequence --> single JSON with all matches for the batch
    Add InterPro signature and entry desc and names
    Add PAINT annotations (if panther is enabled)
    Add go terms (if enabled)
    Add pathways (if enabled)
    Identify representative domains and families
    */
    ch_results = INTERPRO(
        match_results,
        apps,
        data_dir,
        params.xRefsConfig.databases,
        params.xRefsConfig.entries,
        params.xRefsConfig.goterms,
        params.xRefsConfig.pathways,
        params.goterms,
        params.pathways,
        params.appsConfig.panther.paint
    )

    OUTPUT(
        ch_results,
        seq_db_path,
        formats,
        outut_dir,
        params.nucleic,
        workflow.manifest.version
    )
}

workflow.onComplete = {
    def input_file = file(params.input)
    def outputFileName = input_file.getName()
    def outputDir = params.outdir.endsWith('/') ? params.outdir[0..-2] : params.outdir

    println "InterProScan workflow completed successfully: $workflow.success."
    println "Any results are located at ${outputDir}/${outputFileName}.ips6.*"
    println "Duration: $workflow.duration"
}