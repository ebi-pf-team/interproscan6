process PARSER {
    input:
    tuple path(sequences), val(application)
    val preproc_tbl
    val preproc_out

    output:
    tuple path("parse_seq_${sequences}.json"), path("parse_${application}_${sequences}.tbl"), path("parse_${application}_${sequences}.out")

    script:
    """
    python3 $projectDir/scripts/sequences_parse.py -seq ${sequences} > parse_seq_${sequences}.json
    python3 $projectDir/scripts/hmm_parser/hmmer_tbl_parser.py -appl ${application} -preproc_tbl ${preproc_tbl} > parse_${application}_${sequences}.tbl
    python3 $projectDir/scripts/hmm_parser/hmmer_out_parser.py -appl ${application} -preproc ${preproc_out} > parse_${application}_${sequences}.out
    """
}
