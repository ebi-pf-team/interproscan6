nextflow.enable.dsl=2

include { INIT_PIPELINE                 } from "./interproscan/subworkflows/init"
include { PREPARE_SEQUENCES             } from "./interproscan/subworkflows/prepare_sequences"
include { PRECALCULATED_MATCHES         } from "./interproscan/subworkflows/precalculated_matches"
include { SCAN_SEQUENCES                } from "./interproscan/subworkflows/scan"
include { PREPARE_DATA                  } from "./interproscan/subworkflows/prepare_data"
include { OUTPUT                        } from "./interproscan/subworkflows/output"

include { XREFS                         } from "./interproscan/modules/xrefs"
include { REPRESENTATIVE_LOCATIONS      } from "./interproscan/modules/representative_locations"


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
    fasta_file         = Channel.fromPath(INIT_PIPELINE.out.fasta.val)
    outdir             = INIT_PIPELINE.out.outdir.val
    formats            = INIT_PIPELINE.out.formats.val
    applications       = INIT_PIPELINE.out.apps.val
    signalp_mode       = INIT_PIPELINE.out.signalp_mode.val
    interpro_version   = INIT_PIPELINE.out.version.val
    outdir             = INIT_PIPELINE.out.outdir.val
    datadir            = INIT_PIPELINE.out.datadir.val

    PREPARE_DATA(
        applications,
        params.appsConfig,
        interpro_version,
        workflow.manifest.version,
        datadir,
        params.download
    )
    member_db_releases = PREPARE_DATA.out.memberDbReleases.val

    PREPARE_SEQUENCES(fasta_file, applications)
    ch_seqs            = PREPARE_SEQUENCES.out.ch_seqs
    seq_db_path        = PREPARE_SEQUENCES.out.seq_db_path

    matchResults = Channel.empty()
    if (params.offline) {
        SCAN_SEQUENCES(
            ch_seqs,
            member_db_releases,
            applications,
            params.appsConfig,
            datadir
        )
        matchResults = SCAN_SEQUENCES.out
    } else {
        PRECALCULATED_MATCHES(
            ch_seqs,
            applications,
            interpro_version,
            workflow.manifest,
            params.matchesApiUrl,     // from the cmd-offline
            params.lookupService.url, // from confs
        )
        precalculated_matches = PRECALCULATED_MATCHES.out.precalculatedMatches
        no_matches_fastas     = PRECALCULATED_MATCHES.out.noMatchesFasta

        SCAN_SEQUENCES(
            no_matches_fastas,
            member_db_releases,
            applications,
            params.appsConfig,
            datadir
        )

        def expandedScan = SCAN_SEQUENCES.out.flatMap { scan ->
            scan[1].collect { path -> [scan[0], path] }
        }

        combined = precalculated_matches.concat(expandedScan)
        matchResults = combined.groupTuple()
    }
    // matchResults format: [[meta, [member1.json, member2.json, ..., memberN.json]]
    matchResults.view()

    // /* XREFS:
    // Aggregate matches across all members for each sequence --> single JSON with all matches for the batch
    // Add signature and entry desc and names
    // Add PAINT annotations (if panther is enabled)
    // Add go terms (if enabled)
    // Add pathways (if enabled)
    // */
    // XREFS(
    //     matchResults,
    //     apps,
    //     data_dir,
    //     "${data_dir}/interpro/${interpro_release}",
    //     member_db_releases,
    //     params.xRefsConfig.entries,
    //     params.xRefsConfig.goterms,
    //     params.xRefsConfig.pathways,
    //     params.goterms,
    //     params.pathways,
    //     "${data_dir}/${member_db_releases.panther}/${params.appsConfig.panther.paint}"
    // )

    // REPRESENTATIVE_LOCATIONS(XREFS.out)
    // // Collect all JSON files into a single channel so we don't have cocurrent writing to the output files
    // ch_results = REPRESENTATIVE_LOCATIONS.out
    //     .map { meta, json -> json }
    //     .collect()

//     OUTPUT(ch_results)
}

workflow.onComplete = {
    def input_file = file(params.input)
    def outputFileName = input_file.getName()
    def outputDir = params.outdir.endsWith('/') ? params.outdir[0..-2] : params.outdir

    println "InterProScan workflow completed successfully: $workflow.success."
    println "Any results are located at ${outputDir}/${outputFileName}.ips6.*"
    println "Duration: $workflow.duration"
}