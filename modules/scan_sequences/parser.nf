process PARSER {
    input:
    val sequences
    val application
    val preproc

    output:
    path "parse_${application}_${sequences}.out"

    script:
    """
    python $projectDir/scripts/members_parser/${application}.py -seq ${sequences} -preproc ${preproc} > parse_${application}_${sequences}.out
    """
}
