process ESL_TRANSLATE {
    label 'small', 'ips6_container'

    input:
    path fasta

    output:
    path "translated.fasta"

    script:
    """
    /opt/easel/miniapps/esl-translate ${fasta} > translated.fasta
    """
}
