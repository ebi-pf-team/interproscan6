process LOOKUP_CHECK {
    // Checks if protein sequence is included in InterPro
    label 'small'
    
    input:
    val hash_seq
    val is_test

    output:
    path "checked_md5"

    script:
    if ( is_test )
        """
        cat $projectDir/tests/unit_tests/test_outputs/precalc_match_lookup/lookup_check_out > checked_md5
        """
    else
        """
        python3 $projectDir/interproscan/scripts/lookup/lookup_check.py \\
            ${hash_seq} \\
            ${params.url_precalc}${params.check_precalc} \\
            ${params.lookup_retries} \\
            checked_md5
        """
}


process LOOKUP_MATCHES {
    /*
    Retrieves precalculated matches from the Match Lookup Service (MLS).
    A protein sequence can be in the MLS but have no matches associated with it,
    this situation is not detected by LOOKUP_CHECK which only detects that 
    protein sequence has been analysed during an InterPro release.
    */
    label 'small'
    
    input:
    val checked_md5
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
            ${checked_md5} \\
            '${appl}' \\
            ${params.url_precalc}${params.matches} \\
            ${params.lookup_retries} \\
            parsed_match_lookup
        """
}


process LOOKUP_NO_MATCHES {
    label 'small'

    input:
    val checked_md5

    output:
    path "no_match_lookup_fasta.fasta", optional: true

    script:
    """
    output=\$(python3 $projectDir/interproscan/scripts/lookup/lookup_no_matches.py "${checked_md5}")
    
    if [ -n "\$output" ]; then
        echo "\$output" > no_match_lookup_fasta.fasta
    fi
    """
}
