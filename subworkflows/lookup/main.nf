include { LOOKUP_MATCHES } from "../../modules/lookup"

workflow LOOKUP {
    // Prepare connection and retrieve precalculated matched from the InterPro API
    take:
    ch_seqs               // channel of tuples (index, fasta) - fasta files of protein sequences to analyse
    matches_api_apps      // list[str], member db analyses to run that are in the matches API
    url                   // str, url to matches api
    chunk_size            // int
    max_retries           // int

    main:
    LOOKUP_MATCHES(
        ch_seqs.map { index, fasta ->
            tuple(index, fasta, matches_api_apps, url, chunk_size, max_retries)
        }
    )
    precalculatedMatches = LOOKUP_MATCHES.out[0]
    noMatchesFasta = LOOKUP_MATCHES.out[1]

    emit:
    precalculatedMatches
    noMatchesFasta
}
