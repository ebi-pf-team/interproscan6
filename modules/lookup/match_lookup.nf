process MATCHLOOKUP {
    input:
    val fasta
    val applications

    output:
    path 'parsed_match_lookup.out'

    script:
    """
    python $projectDir/scripts/lookup/match_lookup.py -fasta ${fasta} -appl ${applications} > parsed_match_lookup.out
    """
}
