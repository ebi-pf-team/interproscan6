process LOOKUP_MATCHES {
    input:
    val checked_lookup
    val appl

    output:
    path "parsed_match_lookup"

    script:
    """
    python3 $projectDir/scripts/lookup/lookup_matches.py -checked ${checked_lookup} -appl ${appl} -url ${params.url_precalc}${params.matches} > parsed_match_lookup
    """
}
