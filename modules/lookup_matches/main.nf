process LOOKUP_MATCHES {
    input:
    val appl
    val checked_lookup

    output:
    path "parsed_match_lookup_${appl}_${checked_lookup}"

    script:
    """
    python3 $projectDir/scripts/lookup/lookup_matches.py -appl ${appl} -checked ${checked_lookup} -url ${params.lookup.url_precalc}${params.lookup.matches} > parsed_match_lookup_${appl}_${checked_lookup}
    """
}
