process MATCHLOOKUP {
    input:
    val hash_seq
    val appl

    output:
    path 'parsed_match_lookup.out'

    script:
    """
    python3 $projectDir/scripts/lookup/match_lookup.py -seq ${hash_seq} -appl ${appl} > parsed_match_lookup.out
    """
}
