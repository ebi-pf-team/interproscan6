nextflow.enable.dsl=2

include { INIT_PIPELINE                 } from "./interproscan/subworkflows/init"
include { DOWNLOAD_DATA                 } from "./interproscan/subworkflows/download_data"
include { CHECK_DATA                    } from "./interproscan/subworkflows/check_data"
include { PREPARE_SEQUENCES             } from "./interproscan/subworkflows/prepare_sequences"
include { SCAN_SEQUENCES                } from "./interproscan/subworkflows/scan"

include { LOOKUP_MATCHES                } from "./interproscan/modules/lookup"
include { XREFS                         } from "./interproscan/modules/xrefs"
include { REPRESENTATIVE_LOCATIONS      } from "./interproscan/modules/representative_locations"
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
    fasta_file         = Channel.fromPath(INIT_PIPELINE.out.fasta.val)
    outdir             = INIT_PIPELINE.out.outdir.val
    formats            = INIT_PIPELINE.out.formats.val
    apps               = INIT_PIPELINE.out.apps.val
    signalp_mode        = INIT_PIPELINE.out.signalpMode.val
    matches_api_url      = INIT_PIPELINE.out.matchesApiUrl.val

    if (params.download) {
        // Pass the interproRelease to Check data to force sequentialiality 
        DOWNLOAD_DATA(apps)
        downloaded_interpro_release = DOWNLOAD_DATA.out.interproRelease.val

        CHECK_DATA(apps, downloaded_interpro_release)
        data_dir           = CHECK_DATA.out.datadir.val
        interpro_release   = CHECK_DATA.out.interproRelease.val
        member_db_releases = CHECK_DATA.out.memberDbReleases.val
    } else {
        interpro_placeholder = ""
        CHECK_DATA(apps, interpro_placeholder)
        data_dir           = CHECK_DATA.out.datadir.val
        interpro_release   = CHECK_DATA.out.interproRelease.val
        member_db_releases = CHECK_DATA.out.memberDbReleases.val
    }

    PREPARE_SEQUENCES(fasta_file, apps)
    ch_seqs            = PREPARE_SEQUENCES.out.ch_seqs
    seq_db_path        = PREPARE_SEQUENCES.out.seq_db_path

    matchResults = Channel.empty()
    if (matches_api_url != null) {
        LOOKUP_MATCHES(
            ch_seqs,
            apps,
            matches_api_url,
            params.lookupService.chunkSize,
            params.lookupService.maxRetries
        )

        SCAN_SEQUENCES(
            LOOKUP_MATCHES.out[1],
            member_db_releases,
            apps,
            params.appsConfig,
            data_dir
        )

        def expandedScan = SCAN_SEQUENCES.out.flatMap { scan ->
            scan[1].collect { path -> [scan[0], path] }
        }

        def combined = LOOKUP_MATCHES.out[0].concat(expandedScan)
        matchResults = combined.groupTuple()
    } else {
        SCAN_SEQUENCES(
            ch_seqs,
            member_db_releases,
            apps,
            params.appsConfig,
            data_dir
        )
        matchResults = SCAN_SEQUENCES.out
    }
    // matchResults format: [[meta, [member1.json, member2.json, ..., memberN.json]]

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
        "${data_dir}/interpro/${interpro_release}",
        member_db_releases,
        params.xRefsConfig.entries,
        params.xRefsConfig.goterms,
        params.xRefsConfig.pathways,
        params.goterms,
        params.pathways,
        "${data_dir}/${member_db_releases.panther}/${params.appsConfig.panther.paint}"
    )

    REPRESENTATIVE_LOCATIONS(XREFS.out)
    // Collect all JSON files into a single channel so we don't have cocurrent writing to the output files
    ch_results = REPRESENTATIVE_LOCATIONS.out
        .map { meta, json -> json }
        .collect()

    def fileName = params.input.split('/').last()
    def outFileName = "${params.outdir}/${fileName}"

    if (formats.contains("JSON")) {
        WRITE_JSON_OUTPUT(ch_results, "${outFileName}", seq_db_path, params.nucleic, workflow.manifest.version)
    }
    if (formats.contains("TSV")) {
        WRITE_TSV_OUTPUT(ch_results, "${outFileName}", seq_db_path, params.nucleic)
    }
    if (formats.contains("XML")) {
        WRITE_XML_OUTPUT(ch_results, "${outFileName}", seq_db_path, params.nucleic, workflow.manifest.version)
    }
}

workflow.onComplete = {
    def input_file = file(params.input)
    def outputFileName = input_file.getName()
    def outputDir = params.outdir.endsWith('/') ? params.outdir[0..-2] : params.outdir

    println "InterProScan workflow completed successfully: $workflow.success."
    println "Any results are located at ${outputDir}/${outputFileName}.ips6.*"
    println "Duration: $workflow.duration"
}