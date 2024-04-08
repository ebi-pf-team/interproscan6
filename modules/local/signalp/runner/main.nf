process SIGNALP_RUNNER {
    container 'docker.io/library/signalp6'

    input:
    tuple path(fasta), val(mode), path(model_dir), val(organism), val(switches)

    output:
    path "signalp_out"

    script:
    """
    signalp6 --organism ${organism} --fastafile ${fasta} --output_dir signalp_out --mode ${mode}
    """
}
