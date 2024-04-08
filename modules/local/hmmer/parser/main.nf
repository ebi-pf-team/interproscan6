process HMMER_PARSER {
    input:
    path out
    path domtbl
    val tsv_pro

    output:
    path "hmmer_parsed_*"

    script:
    """
    python3 $projectDir/scripts/hmmer/${tsv_pro ? "parser_out" : "parser_domtbl"}.py ${tsv_pro ? "${out}" : "${domtbl}"} > hmmer_parsed_${out}.json
    """
}
