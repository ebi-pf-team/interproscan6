process HMMER_RUNNER {
    container 'docker.io/biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1'

    input:
    tuple path(fasta), path(hmm), val(switches)

    output:
    path "hmmer_${hmm}.out"
    path "hmmer_${hmm}.dtbl"

    script:
    """
    hmmsearch ${switches} -o hmmer_${hmm}.out --domtblout hmmer_${hmm}.dtbl ${hmm} ${fasta}
    """
}
