// codenarc-disable ModuleIncludedTwiceRule
include { PREPARE_APPLICATIONS } from "../subworkflows/prepare/applications"
include { PREPARE_DATABASES    } from "../subworkflows/prepare/databases"
include { PREPARE_SEQUENCES    } from "../subworkflows/prepare/sequences"
include { LOOKUP               } from "../subworkflows/lookup"
include { SCAN_SEQUENCES as SCAN_REMAINING;
          SCAN_SEQUENCES as SCAN_LOCALLY;
          SCAN_SEQUENCES       } from "../subworkflows/scan"
include { COMBINE              } from "../subworkflows/combine"
include { OUTPUT               } from "../subworkflows/output"

import uk.ac.ebi.interpro.InterProScan

workflow INTERPROSCAN {
    take:
    fasta_file            // Channel.fromPath(input fasta file)
    applications          // list[str], names of applications to run
    apps_config           // map, contents of the conf/applications.conf file
    data_dir              // str, path to the data directory
    outprefix             // str, prefix for output files
    formats               // list[str], output formats
    interpro_version      // str, version of InterPro
    interproscan_version  // str, version of InterProScan
    interproscan_name     // str, name of the InterProScan workflow: "InterProScan6"
    no_matches_api        // boolean, use the Matches API
    matches_api_url       // str, url to the Matches API
    matches_api_chunk_size  // int, chunk size for Matches API requests
    matches_api_max_retries // int, max retries for Matches API requests
    batch_size            // int, number of sequences to process per batch
    sub_batch_size        // int, number of sequences per subbatch
    nucleic               // boolean, input is nucleic acid sequences
    skip_repr_locations   // boolean, skip indentifying representative locations
    goterms               // boolean, include GO terms in the output files
    pathways              // boolean, include pathway terms in the output files
    globus                // boolean, use Globus to download databases
    enforce_compatibility // boolean, whether to enforce compatibility between IPRScan and InterPro data versions

    main:

    PREPARE_APPLICATIONS(
        applications,
        no_matches_api,
        matches_api_url,
        interproscan_name,
        interproscan_version
    )
    local_only_apps   = PREPARE_APPLICATIONS.out.local_only_apps.val
    matches_api_apps  = PREPARE_APPLICATIONS.out.matches_api_apps.val
    api_version       = PREPARE_APPLICATIONS.out.api_version.val

    PREPARE_DATABASES(
        local_only_apps,
        matches_api_apps,
        apps_config,
        data_dir,
        interpro_version,
        interproscan_version,
        api_version,
        goterms,
        pathways,
        globus,
        enforce_compatibility
    )
    use_matches_api = PREPARE_DATABASES.out.use_matches_api.val
    db_releases = PREPARE_DATABASES.out.versions

    PREPARE_SEQUENCES(
        fasta_file,
        nucleic,
        batch_size
    )
    ch_seqs              = PREPARE_SEQUENCES.out.ch_seqs
    seq_db_path          = PREPARE_SEQUENCES.out.seq_db_path

    match_results = Channel.empty()

    if (!use_matches_api || matches_api_apps.isEmpty()) {
        SCAN_SEQUENCES(
            ch_seqs,
            db_releases,
            applications,
            apps_config,
            data_dir,
            local_only_apps,
            sub_batch_size
        )
        match_results = SCAN_SEQUENCES.out
    } else {
        /* Retrieve precalculated matches from the Match lookup API
        Then run analyses on sequences not listed in the API */
        LOOKUP(
            ch_seqs,
            matches_api_apps,
            db_releases,
            api_version,
            matches_api_url,
            matches_api_chunk_size,
            matches_api_max_retries
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
            sub_batch_size
        )

        SCAN_LOCALLY(
            ch_seqs,
            db_releases,
            local_only_apps,
            apps_config,
            data_dir,
            matches_api_apps + local_only_apps,
            sub_batch_size
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
        goterms,
        pathways,
        apps_config.panther.paint,
        skip_repr_locations,
        batch_size
    )

    output_files = OUTPUT(
        ch_results,
        seq_db_path,
        formats,
        outprefix,
        nucleic,
        interproscan_version,
        db_releases,
        batch_size
    )

    emit:
    output_files
}

workflow PREPARE_INTERPROSCAN {
    // Parses an interproscan applications.config file and prepares it for use
    take:
    apps_config  // str repr of path to the applications.config file
    applications // list[str], names of applications to prepare
    use_gpu      // boolean, use GPU acceleration where applicable

    main:
    (apps_config, warn) = InterProScan.parseAppsConfig(use_gpu, applications, apps_config)
    if (warn) {
        log.warn warn
    }

    emit:
    apps_config
}
