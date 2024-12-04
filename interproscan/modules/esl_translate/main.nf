process ESL_TRANSLATE {
    label 'small'

    input:
    path fasta

    output:
    tuple path("translated.fasta"), path(fasta)

    script:
    """
    /opt/easel/miniapps/esl-translate ${fasta} > translated.fasta
    """
}
