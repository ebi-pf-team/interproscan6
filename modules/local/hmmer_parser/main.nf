process HMMER_PARSER {
    input:
    path hmmer_result

    output:
    path "hmmer_parsed.json"

    script:
    """
    python3 $projectDir/scripts/hmm_parser/hmmer_parser_domtbl.py hmmer_result > hmmer_parsed.json
    """
}
