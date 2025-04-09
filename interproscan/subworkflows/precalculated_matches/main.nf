include { LOOKUP_MATCHES;
          getMatchesApiUrl } from "../../modules/lookup"

workflow PRECALCULATED_MATCHES {
    // Prepare connection and retrive precalculated matched from the InterPro API
    take:
    ch_seqs           // fasta files of protein sequences to analyse
    apps              // member db analyses to run
    member_db_releases   // str, interpro db version number
    iprscanRelease    // str, full iprscan release number
    workflowManifest  // map
    matchesApiUrl     // str, from cmd-line
    lookupService     // Map of confs/lookup.conf

    main:
    def _matchesApiUrl = getMatchesApiUrl(
        matchesApiUrl, lookupService.url, member_db_releases.interpro, iprscanRelease, workflowManifest, log
    )

    LOOKUP_MATCHES(
        ch_seqs,
        apps,
        _matchesApiUrl,
        lookupService.chunkSize,
        lookupService.maxRetries
    )

    precalculatedMatches = LOOKUP_MATCHES.out[0]
    noMatchesFasta       = LOOKUP_MATCHES.out[1]

    emit:
    precalculatedMatches
    noMatchesFasta

}