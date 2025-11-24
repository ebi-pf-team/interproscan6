process ESL_TRANSLATE {
    label 'mem_low', 'time_short', 'ips6_container'

    input:
    path fasta

    output:
    path "translated.fasta"

    script:
    """
    esl-translate ${fasta} > translated.fasta
    """
}
