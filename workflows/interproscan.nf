include { PREPARE_APPLICATIONS } from "${moduleDir}/../subworkflows/prepare/applications"
include { PREPARE_DATABASES    } from "${moduleDir}/../subworkflows/prepare/databases"
include { PREPARE_SEQUENCES    } from "${moduleDir}/../subworkflows/prepare/sequences"
include { LOOKUP               } from "${moduleDir}/../subworkflows/lookup"
include { SCAN_SEQUENCES as SCAN_REMAINING;
          SCAN_SEQUENCES as SCAN_LOCALLY;
          SCAN_SEQUENCES       } from "${moduleDir}/../subworkflows/scan"
include { COMBINE              } from "${moduleDir}/../subworkflows/combine"
include { OUTPUT               } from "${moduleDir}/../subworkflows/output"

workflow INTERPROSCAN {
    take:
    fasta_file            // Channel.fromPath(input fasta file)
    applications          // list[str], names of applications to run
    no_matches_api        // boolean, whether to use the Matches API
    matches_api_url       // str, url to the Matches API
    apps_config           // map, contents of the conf/applications.conf file
    data_dir              // str, path to the data directory
    outprefix             // str, prefix for output files
    formats               // list[str], output formats
    interpro_version      // str, version of InterPro
    interproscan_version  // str, version of InterProScan
    interproscan_name     // str, name of the InterProScan workflow: "InterProScan6"

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
        interproscan_version,
        db_releases,
        params.batchSize
    )

}