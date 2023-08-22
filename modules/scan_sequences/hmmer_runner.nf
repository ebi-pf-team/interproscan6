process HMMER_RUNNER {
    input:
    val fasta_path
    val hmm_path
    val options

    script:
    """
    echo ${params.hmmsearch_bin} ${hmm_path} ${fasta_path} ${options} ${hmm_path}_${fasta_path}.out
    """
}

//     output:
//     path "${hmm_path}_${fasta_path}.out"
// ${params.hmmsearch_bin} ${hmm_path} ${fasta_path} ${options} > ${hmm_path}_${fasta_path}.out