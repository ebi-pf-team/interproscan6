process HMMER_RUNNER {
    container 'interproscan6'

    input:
    tuple val(fasta), path(hmm), val(switches)

    output:
    path "hmmer_${hmm}.out"
    path "hmmer_${hmm}.dtbl"

    script:
    """
    ${params.hmmsearch_bin} ${switches} -o hmmer_${hmm}.out --domtblout hmmer_${hmm}.dtbl ${hmm} ${fasta}
    """
}
