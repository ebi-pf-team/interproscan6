process HMMER_PARSER {
    input:
    path preproc_out
    path preproc_domtbl
    val tsv_pro

    output:
    path "hmm_parsed_${preproc_out}.json"

    script:
    """
    python3 $projectDir/scripts/hmm_parser/${tsv_pro ? "hmmer_parser_out" : "hmmer_parser_domtbl"}.py -hmmer_file ${tsv_pro ? "${preproc_out}" : "${preproc_domtbl}"} > hmm_parsed_${preproc_out}.json
    """
}
