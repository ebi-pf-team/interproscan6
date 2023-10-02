process PARSER {
    input:
    path preproc_out
    path preproc_domtbl

    output:
    path "hmm_parsed_${preproc_domtbl}.json"

    script:
    """
    python3 $projectDir/scripts/hmm_parser/hmmer_parser.py -out ${preproc_out} -domtbl ${preproc_domtbl} > hmm_parsed_${preproc_domtbl}.json
    """
}
