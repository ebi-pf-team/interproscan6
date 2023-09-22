process PARSER {
    input:
    val preproc_out
    val preproc_domtbl

    output:
    path "hmm_parsed.json"

    script:
    """
    python3 $projectDir/scripts/hmm_parser/hmmer_parser.py -out ${preproc_out} -domtbl ${preproc_domtbl} > hmm_parsed.json
    """
}
