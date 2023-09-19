process PARSER {
    input:
    tuple path(sequences), val(application)
    val preproc_out
    val preproc_tbl

    output:
    path "parsed_result_${application}.json"

    script:
    """
    python3 $projectDir/scripts/sequences_parse.py -seq ${sequences} > parsed_seq.json
    python3 $projectDir/scripts/hmm_parser/hmmer_tbl_parser.py -appl ${application} -preproc ${preproc_tbl} > parsed_${application}.tbl
    python3 $projectDir/scripts/hmm_parser/hmmer_out_parser.py -appl ${application} -preproc ${preproc_out} > parsed_${application}.out
    python3 $projectDir/scripts/hmm_parser/build_parsed_output.py -seq parsed_seq.json -tbl parsed_${application}.tbl -out parsed_${application}.out > parsed_result_${application}.json
    """
}
