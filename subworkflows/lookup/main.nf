include { check_matches_api; GET_MATCHES } from "../../modules/lookup"

workflow LOOKUP {
    take:
    fasta                 // [meta, fasta]
    applications          // list[str]
    matches_api_url       // str, url to matches api
    chunk_size            // int
    max_retries           // int
    interproscan_version  // str

    main:
    (appls_in_api, appls_not_in_api) = check_matches_api(applications, matches_api_url, interproscan_version)

    GET_MATCHES(
        fasta,
        appls_in_api,
        matches_api_url,
        chunk_size,
        max_retries
    )

    emit:
    json             = GET_MATCHES.out.json
    fasta            = GET_MATCHES.out.fasta
    appls_in_api
    appls_not_in_api
}
