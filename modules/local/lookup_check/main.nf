process LOOKUP_CHECK {
    input:
    val hash_seq

    output:
    path checked_precalc_matches

    script:
    """
    python3 $projectDir/scripts/lookup/lookup_check.py -seq ${hash_seq} -url ${params.url_precalc}${params.check_precalc} > checked_precalc_matches
    """
}
