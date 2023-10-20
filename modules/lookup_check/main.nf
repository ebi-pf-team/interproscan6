process LOOKUP_CHECK {
    input:
    val hash_seq

    output:
    path "no_match_${hash_seq}"

    script:
    """
    python3 $projectDir/scripts/lookup/lookup_check.py -seq ${hash_seq} -url ${params.lookup.url_precalc}${params.lookup.check_precalc} > no_match_${hash_seq}
    """
}
