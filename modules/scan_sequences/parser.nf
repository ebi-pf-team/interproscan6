process PARSER {
    input:
    tuple path(sequences), val(application)
    val preproc

    output:
    path "parse_${application}_${sequences}.out"

    script:
    """
    python3 $projectDir/scripts/hmm_parser/hmmer_parser.py -appl ${application} -seq ${sequences} -preproc ${preproc} > parse_${application}_${sequences}.out
    """
}
