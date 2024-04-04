process HMMER_RUNNER {
    container 'docker.io/staphb/hmmer:latest'

    input:
    tuple path(fasta), path(hmm), val(switches)

    output:
    path "hmmer_${hmm}.out"
    path "hmmer_${hmm}.dtbl"

    script:
    """
    ${params.bin.hmmer.hmmsearch} ${switches} -o hmmer_${hmm}.out --domtblout hmmer_${hmm}.dtbl ${hmm} ${fasta}
    """
}
