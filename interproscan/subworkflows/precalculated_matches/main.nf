include { LOOKUP_MATCHES;
          getMatchesApiUrl } from "../../modules/lookup"

workflow PRECALCULATED_MATCHES {
    // Prepare connection and retrive precalculated matched from the InterPro API
    take:
    ch_seqs           // fasta files of protein sequences to analyse
    apps              // member db analyses to run
    interproRelease   // str, interpro db version number
    iprscanRelease    // str, full iprscan release number
    matchesApiUrl     // str, from cmd-line
    lookupServiceUrl  // str, from confs/lookup.conf

    main:
    _matchesApiUrl = getMatchesApiUrl(
        matchesApiUrl, lookupServiceUrl, interproRelease, iprscanRelease, log
    )

    LOOKUP_MATCHES(
        ch_seqs,
        apps,
        _matchesApiUrl,
        params.lookupService.chunkSize,
        params.lookupService.maxRetries
    )

    precalculatedMatches = LOOKUP_MATCHES.out[0]
    noMatchesFasta       = LOOKUP_MATCHES.out[1]

    emit:
    precalculatedMatches
    noMatchesFasta

}