process HMMER_PARSER {
    input:
    path preproc_out
    path preproc_domtbl
    val tsv_pro

    output:
    path "hmm_parsed_${preproc}.json"

    script:
    """
    def hmmer_file = tsv_pro ? ${preproc_out} : ${preproc_domtbl}
    python3 $projectDir/scripts/hmm_parser/hmmer_parser.py -hmmer_file ${hmmer_file} > hmm_parsed_${preproc}.json
    """
}
