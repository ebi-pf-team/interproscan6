process ESL_TRANSLATE {
    label     'mem_low'
    label     'time_short'
    container 'interpro/esl-translate:0.46'

    input:
    path fasta

    output:
    path "translated.fasta"

    script:
    """
    esl-translate ${fasta} > translated.fasta
    """
}
