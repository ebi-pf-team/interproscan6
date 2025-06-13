nextflow.enable.dsl=2

include { INIT_PIPELINE      } from "./subworkflows/init"
include { PREPARE_DATABASES  } from "./subworkflows/prepare/databases"
include { PREPARE_SEQUENCES  } from "./subworkflows/prepare/sequences"
include { LOOKUP             } from "./subworkflows/lookup"
include { SCAN_SEQUENCES     } from "./subworkflows/scan"
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
        params.datadir,
        params.formats,
        params.outdir,
        params.outprefix,
        params.noMatchesApi,
        params.interpro,
        params.skipInterpro,
        params.skipApplications,
        params.goterms,
        params.pathways
    )
    fasta_file           = Channel.fromPath(INIT_PIPELINE.out.fasta.val)
    applications         = INIT_PIPELINE.out.apps.val
    data_dir             = INIT_PIPELINE.out.datadir.val
    outprefix            = INIT_PIPELINE.out.outprefix.val
    formats              = INIT_PIPELINE.out.formats.val
    interpro_version     = INIT_PIPELINE.out.version.val

    PREPARE_DATABASES(
        applications,
        params.appsConfig,
        data_dir,
        interpro_version,
        workflow.manifest.version,
        params.goterms,
        params.pathways
    )
    db_releases = PREPARE_DATABASES.out.versions
    interproscan_version = PREPARE_DATABASES.out.iprscan_major_minor

    PREPARE_SEQUENCES(
        fasta_file,
        applications
    )
    ch_seqs              = PREPARE_SEQUENCES.out.ch_seqs
    seq_db_path          = PREPARE_SEQUENCES.out.seq_db_path

    match_results = Channel.empty()

    if (params.noMatchesApi) {
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
        LOOKUP(
            ch_seqs,
            applications,
            db_releases,
            interproscan_version,
            workflow.manifest,
            params.matchesApiUrl,
            params.matchesApiChunkSize,
            params.matchesApiMaxRetries
        )
        precalculated_matches = LOOKUP.out.precalculatedMatches
        no_matches_fastas     = LOOKUP.out.noMatchesFasta

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
        params.appsConfig.panther.paint,
        params.skipInterpro
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
}
