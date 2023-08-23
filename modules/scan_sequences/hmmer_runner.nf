process HMMER_RUNNER {
    input:
    tuple path(sequences), path(hmm), val(switches)

    output:
    path "${hmm_path}_${fasta_path}.out"

    script:
    """
    ${params.hmmsearch_bin} ${hmm_path} ${fasta_path} ${options} > ${hmm_path}_${fasta_path}.out
    """
}
