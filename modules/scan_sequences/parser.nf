process PARSER {
    input:
    tuple path(sequences), val(application)
    val preproc

    output:
    path "parse_${application}_${sequences}.out"

    script:
    """
    python $projectDir/scripts/members_parser/hmmer_parser.py -appl ${application} -seq ${sequences} -preproc ${preproc} > parse_${application}_${sequences}.out
    """
}
