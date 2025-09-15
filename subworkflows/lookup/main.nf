include { PREPARE_LOOKUP; LOOKUP_MATCHES } from "../../modules/lookup"

workflow LOOKUP {
    // Prepare connection and retrieve precalculated matched from the InterPro API
    take:
    ch_seqs               // fasta files of protein sequences to analyse
    matches_api_apps      // member db analyses to run that are in the matches API
    db_releases           // map: [db: version, dirpath]           
    interproscan_version  // major.minor interproscan version number
    api_version           // version of the matches API
    url                   // str, url to matches api
    chunk_size            // int
    max_retries           // int

    main:
    PREPARE_LOOKUP(
        matches_api_apps,
        api_version,
        db_releases,
        url
    )

    // Branch sequences based on API availability
    api_result = PREPARE_LOOKUP.out[0]
        .combine(ch_seqs)
        .branch {
            available: it[0] != null
            unavailable: it[0] == null
        }

    // Run LOOKUP_MATCHES only on available branch
    LOOKUP_MATCHES(
        api_result.available.map { api_url, index, fasta ->
            tuple(index, fasta, matches_api_apps, api_url, chunk_size, max_retries)
        }
    )

    precalculatedMatches = LOOKUP_MATCHES.out[0]
        .mix(
            api_result.unavailable.map { _, index, fasta -> 
                tuple(index, null)
            }
        )

    noMatchesFasta = LOOKUP_MATCHES.out[1]
        .mix(
            api_result.unavailable.map { _, index, fasta -> 
                tuple(index, fasta)
            }
        )

    emit:
    precalculatedMatches
    noMatchesFasta
}
