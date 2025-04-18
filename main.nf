nextflow.enable.dsl=2
import groovy.json.JsonSlurper  // until selective downloads

include { INIT_PIPELINE                 } from "./interproscan/subworkflows/init"
include { PREPARE_DATA                  } from "./interproscan/subworkflows/prepare_data"
include { PREPARE_SEQUENCES             } from "./interproscan/subworkflows/prepare_sequences"
include { PRECALCULATED_MATCHES         } from "./interproscan/subworkflows/precalculated_matches"
include { SCAN_SEQUENCES                } from "./interproscan/subworkflows/scan"
include { INTERPRO                      } from "./interproscan/subworkflows/interpro"
include { OUTPUT                        } from "./interproscan/subworkflows/output"

workflow {
    println "# ${workflow.manifest.name} ${workflow.manifest.version}"
    println "# ${workflow.manifest.description}\n"

    if (params.keySet().any { it.equalsIgnoreCase("help") }) {
        InterProScan.printHelp(params.appsConfig)
        exit 0
    }

    // Params validation
    InterProScan.validateParams(params, log)

    INIT_PIPELINE(
        params.input,
        params.applications,
        params.appsConfig,
        params.download,
        params.offline,
        params.datadir,
        params.formats,
        params.outdir,
        params.signalpMode,
        params.matchesApiUrl,
        params.interpro
    )
    fasta_file           = Channel.fromPath(INIT_PIPELINE.out.fasta.val)
    applications         = INIT_PIPELINE.out.apps.val
    data_dir             = INIT_PIPELINE.out.datadir.val
    out_dir              = INIT_PIPELINE.out.outdir.val
    formats              = INIT_PIPELINE.out.formats.val
    signalp_mode         = INIT_PIPELINE.out.signalp_mode.val
    interpro_version     = INIT_PIPELINE.out.version.val

    PREPARE_DATA(
        applications,
        params.appsConfig,
        data_dir,
        interpro_version,
        workflow.manifest.version,
        params.download,
        params.goterms,
        params.pathways
    )
    db_releases   = PREPARE_DATA.out.versions
    interproscan_version = PREPARE_DATA.out.iprscan_major_minor

    PREPARE_SEQUENCES(
        fasta_file,
        applications
    )
    ch_seqs              = PREPARE_SEQUENCES.out.ch_seqs
    seq_db_path          = PREPARE_SEQUENCES.out.seq_db_path

    match_results = Channel.empty()

    if (params.offline) {
        SCAN_SEQUENCES(
            ch_seqs,
            db_releases,
            applications,
            params.appsConfig,
            data_dir
        )
        match_results = SCAN_SEQUENCES.out
    } else {
        /* Retrieve precalculated matches from the Match lookup API
        Then run analyses on sequences not listed in the API */
        PRECALCULATED_MATCHES(
            ch_seqs,
            applications,
            db_releases,
            interproscan_version,
            workflow.manifest,
            params.matchesApiUrl,     // from the cmd-offline
            params.lookupService,     // from confs
        )
        precalculated_matches = PRECALCULATED_MATCHES.out.precalculatedMatches
        no_matches_fastas     = PRECALCULATED_MATCHES.out.noMatchesFasta

        SCAN_SEQUENCES(
            no_matches_fastas,
            db_releases,
            applications,
            params.appsConfig,
            data_dir
        )

        def expandedScan = SCAN_SEQUENCES.out.flatMap { scan ->
            scan[1].collect { path -> [scan[0], path] }
        }

        combined = precalculated_matches.concat(expandedScan)
        match_results = combined.groupTuple()
    }
    // match_results format: [[meta, [member1.json, member2.json, ..., memberN.json]]

    /* INTERPRO:
    Aggregate matches across all members for each sequence --> single JSON with all matches for the batch
    Add InterPro signature and entry desc and names, PAINT annotations (panther only),
    go terms (if enabled), and pathways (if enabled). Then identify representative domains and families
    */
    ch_results = INTERPRO(
        match_results,
        applications,
        db_releases,
        params.goterms,
        params.pathways,
        params.appsConfig.panther.paint
    )

    OUTPUT(
        ch_results,
        seq_db_path,
        formats,
        out_dir,
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