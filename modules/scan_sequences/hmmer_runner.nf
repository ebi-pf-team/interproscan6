process HMMER_RUNNER {
    input:
    tuple path(sequences), path(hmm), val(switches)

    output:
    path "${hmm_path}_${fasta_path}.out"

    script:
    """
    ${params.hmmsearch_bin} ${hmm} ${sequences} ${switches} > ${hmm}_${sequences}.out
    """
}
