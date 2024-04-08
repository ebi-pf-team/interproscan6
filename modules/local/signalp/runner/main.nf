process SIGNALP_RUNNER {
    input:
    tuple path(fasta), val(mode), path(model_dir), val(organism), val(switches)

    output:
    path "signalp_out"

    script:
    """
    signalp6 --organism ${organism} --fastafile ${fasta} --output_dir signalp_out ${switches} --mode ${mode} --model_dir ${model_dir}
    """
}
