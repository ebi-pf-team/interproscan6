nextflow.enable.dsl=2
import java.time.format.DateTimeFormatter

include { INIT_PIPELINE      } from "./subworkflows/init"
include { PREPARE_DATABASES  } from "./subworkflows/prepare/databases"
include { PREPARE_SEQUENCES  } from "./subworkflows/prepare/sequences"
include { LOOKUP             } from "./subworkflows/lookup"
include { SCAN_SEQUENCES as SCAN_REMAINING;
          SCAN_SEQUENCES as SCAN_LOCALLY;
          SCAN_SEQUENCES     } from "./subworkflows/scan"
include { COMBINE            } from "./subworkflows/combine"
include { OUTPUT             } from "./subworkflows/output"

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
        params.runMl,
        params.useGpu,
        params.datadir,
        params.formats,
        params.outdir,
        params.outprefix,
        params.noMatchesApi,
        params.matchesApiUrl,
        params.interpro,
        params.skipInterpro,
        params.skipApplications,
        params.goterms,
        params.pathways,
        workflow.manifest
    )
    fasta_file           = Channel.fromPath(INIT_PIPELINE.out.fasta.val)
    apps_config          = INIT_PIPELINE.out.apps_config.val
    local_only_apps      = INIT_PIPELINE.out.local_only_apps.val
    matches_api_apps     = INIT_PIPELINE.out.matches_api_apps.val
    api_version          = INIT_PIPELINE.out.api_version.val
    data_dir             = INIT_PIPELINE.out.datadir.val
    outprefix            = INIT_PIPELINE.out.outprefix.val
    formats              = INIT_PIPELINE.out.formats.val
    interpro_version     = INIT_PIPELINE.out.version.val

    PREPARE_DATABASES(
        local_only_apps,
        matches_api_apps,
        apps_config,
        data_dir,
        interpro_version,
        workflow.manifest.version,
        params.noMatchesApi,
        params.goterms,
        params.pathways,
        params.globus
    )
    db_releases = PREPARE_DATABASES.out.versions
    interproscan_version = PREPARE_DATABASES.out.iprscan_major_minor

    PREPARE_SEQUENCES(
        fasta_file,
        params.nucleic,
        params.batchSize
    )
    ch_seqs              = PREPARE_SEQUENCES.out.ch_seqs
    seq_db_path          = PREPARE_SEQUENCES.out.seq_db_path

    match_results = Channel.empty()

    if (params.noMatchesApi || matches_api_apps.isEmpty()) {
        SCAN_SEQUENCES(
            ch_seqs,
            db_releases,
            local_only_apps,
            apps_config,
            data_dir,
            local_only_apps,
            params.subBatchSize
        )
        match_results = SCAN_SEQUENCES.out
    } else {
        /* Retrieve precalculated matches from the Match lookup API
        Then run analyses on sequences not listed in the API */
        LOOKUP(
            ch_seqs,
            matches_api_apps,
            db_releases,
            interproscan_version,
            api_version,
            params.matchesApiUrl,
            params.matchesApiChunkSize,
            params.matchesApiMaxRetries
        )
        precalculated_matches = LOOKUP.out.precalculatedMatches
        no_matches_fastas     = LOOKUP.out.noMatchesFasta

        SCAN_REMAINING(
            no_matches_fastas,
            db_releases,
            matches_api_apps,
            apps_config,
            data_dir,
            matches_api_apps + local_only_apps,
            params.subBatchSize
        )

        SCAN_LOCALLY(
            ch_seqs,
            db_releases,
            local_only_apps,
            apps_config,
            data_dir,
            matches_api_apps + local_only_apps,
            params.subBatchSize
        )

        def expandedRemainingScan = SCAN_REMAINING.out.flatMap { scan ->
            scan[1].collect { path -> [scan[0], path] }
        }
        def expandedLocalScan     = SCAN_LOCALLY.out.flatMap { scan ->
            scan[1].collect { path -> [scan[0], path] }
        }

        def allExpandedScans = expandedRemainingScan.concat(expandedLocalScan)
        combined             = precalculated_matches.concat(allExpandedScans)
        match_results        = combined.groupTuple()
    }
    // match_results format: [[meta, [member1.json, member2.json, ..., memberN.json]]

    /* COMBINE:
    Aggregate matches across all members for each sequence --> single JSON with all matches for the batch
    Add InterPro signature and entry desc and names, PAINT annotations (panther only),
    go terms (if enabled), and pathways (if enabled). Then identify representative domains and families
    */
    ch_results = COMBINE(
        match_results,
        db_releases,
        params.goterms,
        params.pathways,
        apps_config.panther.paint,
        params.skipInterpro,
        params.batchSize
    )

    OUTPUT(
        ch_results,
        seq_db_path,
        formats,
        outprefix,
        params.nucleic,
        workflow.manifest.version,
        db_releases
    )

    workflow.onComplete = {
        if (workflow.success) {
            println "\nInterProScan completed successfully"
            println "Results available at: ${outprefix}.*"
            if (workflow.duration.toSeconds() <= 60) {
                DateTimeFormatter formatter = DateTimeFormatter.ofPattern("dd-MMM-yyyy HH:mm:ss");
                println "Completed at        : ${workflow.complete.format(formatter)}"
                println "Duration            : ${workflow.duration}"
            }
        }        
    }
}
