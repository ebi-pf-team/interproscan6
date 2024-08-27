process LOOKUP_MATCHES {
    label 'io'

    input:
    val checked_lookup
    val appl
    val is_test

    output:
    path("parsed_match_lookup"), optional: true

    /*
    'parsed_match_lookup' will only be generated if any matches are returned from
    the MLS. If no matches are retrieved from the MLS, 'parsed_match_lookup'
    will not be created.
    */
    script:
    if ( is_test )
        """
        cat $projectDir/tests/unit_tests/test_outputs/precalc_match_lookup/lookup_matches_out > parsed_match_lookup
        """
    else
        """
        python3 $projectDir/interproscan/scripts/lookup/lookup_matches.py \\
            ${checked_lookup} \\
            '${appl}' \\
            ${params.url_precalc}${params.matches} \\
            ${params.lookup_retries} \\
            parsed_match_lookup
        """
}
