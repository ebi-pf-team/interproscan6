process HMMER_RUNNER {
    container 'docker.io/biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1'

    input:
    tuple path(fasta), path(hmm), val(switches)

    output:
    path "${hmm}.out"
    path "${hmm}.dtbl"

    script:
    """
    hmmsearch ${switches} -o ${hmm}.out --domtblout ${hmm}.dtbl ${hmm} ${fasta}
    """
}
