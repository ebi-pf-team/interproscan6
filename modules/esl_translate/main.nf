process ESL_TRANSLATE {
    label     'mem_low'
    label     'time_short'
    container 'interpro/esl-translate:0.46'

    input:
    path fasta

    output:
    path "translated.fasta", emit: results
    tuple val(task.process), val('easel'), eval("software-version '^# Easel ([^ ]+).*' esl-translate -h"), topic: versions

    script:
    """
    esl-translate ${fasta} > translated.fasta
    """
}
