process PARSER {
    input:
    val fasta
    val application
    val preproc_output

    output:
    path parser_result

    script:
    """
    echo parse $application ${fasta} ${preproc_output} > parser_result
    """
}

// python $projectDir/scripts/members_parser/${application}.py -fasta ${fasta} -preproc ${preproc_output} > parser_result
