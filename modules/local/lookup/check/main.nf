process LOOKUP_CHECK {
    input:
    val hash_seq

    output:
    path checked_md5

    script:
    """
    python3 $projectDir/scripts/lookup/lookup_check.py ${hash_seq} ${params.url_precalc}${params.check_precalc} > checked_md5
    """
}
