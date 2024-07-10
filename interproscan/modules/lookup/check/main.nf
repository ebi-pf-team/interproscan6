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
        python3 $projectDir/interproscan/scripts/lookup/lookup_check.py ${hash_seq} ${params.url_precalc}${params.check_precalc} > checked_md5
        """
}
