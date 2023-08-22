process PARSER {
    input:
    val fasta
    val application
    val preproc_output

    script:
    """
    echo parse $application ${fasta} ${preproc_output}
    """
}

// python $projectDir/scripts/members_parser/${application}.py -fasta ${fasta} -preproc ${preproc_output} > parsed_hmm
