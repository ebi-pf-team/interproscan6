process HMMER_RUNNER {
    container 'interproscan6'

    input:
    tuple path(sequences), path(hmm), val(switches)

    output:
    path "hmmer_${hmm}_${sequences}.out"

    script:
    """
    ${params.hmmsearch_bin} ${switches} ${hmm} ${sequences} > hmmer_${hmm}_${sequences}.out
    """
}
