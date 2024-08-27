process LOOKUP_CHECK {
    label 'io'

    input:
    val hash_seq
    val is_test

    output:
    path checked_md5

    script:
    if ( is_test )
        """
        cat $projectDir/tests/unit_tests/test_outputs/precalc_match_lookup/lookup_check_out > checked_md5
        """
    else
        """
        python3 $projectDir/interproscan/scripts/lookup/lookup_check.py ${hash_seq} \\
        ${params.url_precalc}${params.check_precalc} ${params.lookup_retries} > checked_md5
        """
}


process LOOKUP_MATCHES {
    label 'io'

    input:
    val checked_md5
    val appl
    val is_test

    output:
    path("parsed_match_lookup"), optional: true

    script:
    if ( is_test )
        """
        cat $projectDir/tests/unit_tests/test_outputs/precalc_match_lookup/lookup_matches_out > parsed_match_lookup
        """
    else
        """
        python3 $projectDir/interproscan/scripts/lookup/lookup_matches.py ${checked_md5} '${appl}' \\
        ${params.url_precalc}${params.matches} ${params.lookup_retries} > result

        if [[ -s result ]]; then
            mv result parsed_match_lookup
        fi
        """
}


process LOOKUP_NO_MATCHES {
    label 'io'

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
