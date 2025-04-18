include { PREPARE_LOOKUP; LOOKUP_MATCHES } from "../../modules/lookup"

workflow PRECALCULATED_MATCHES {
    // Prepare connection and retrieve precalculated matched from the InterPro API
    take:
    ch_seqs               // fasta files of protein sequences to analyse
    apps                  // member db analyses to run
    db_releases           // map: [db: version, dirpath]           
    interproscan_version  // major.minor interproscan version number
    workflow_manifest     // map, from nextflow.conf
    matches_api_url       // str, from cmd-line
    lookup_service        // Map of confs/lookup.conf

    main:
    // Initialise channels for outputs
    precalculatedMatches = Channel.empty()
    noMatchesFasta = Channel.empty()

    _url = matches_api_url ?: lookup_service.url
    PREPARE_LOOKUP(
        _url,
        db_releases,
        interproscan_version,
        workflow_manifest
    )
    _matches_api_url = PREPARE_LOOKUP.out[0]

    _matches_api_url
        .filter { it }
        .combine(ch_seqs)
        .map { url, index, fasta ->
            tuple(index, fasta, apps, url, lookup_service.chunkSize, lookup_service.maxRetries)
        }
        .set { lookup_input }

    LOOKUP_MATCHES(lookup_input)

    precalculatedMatches = LOOKUP_MATCHES.out[0]
    noMatchesFasta = LOOKUP_MATCHES.out[1]

    // fallback when no API URL is available
    _matches_api_url
        .filter { !it }
        .map { _ ->
            precalculatedMatches = Channel.empty()
            noMatchesFasta = ch_seqs
        }

    emit:
    precalculatedMatches
    noMatchesFasta
}