process LOOKUP_MATCHES {
    label 'io'
    
    input:
    val checked_lookup
    val appl
    val is_test

    output:
    path "parsed_match_lookup"

    script:
    if ( is_test )
        """
        echo '\$(cat tests/unit_tests/test_outputs/precalc_match_lookup/lookup_matches_out)' > parsed_match_lookup
        """
    else
        """
        python3 $projectDir/scripts/lookup/lookup_matches.py ${checked_lookup} '${appl}' ${params.url_precalc}${params.matches} > parsed_match_lookup
        """
}
