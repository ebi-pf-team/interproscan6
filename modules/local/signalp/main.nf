process SIGNALP_RUNNER {
    container 'docker.io/library/signalp6'

    input:
    tuple path(fasta), val(mode), path(model_dir), val(organism), val(switches), val(pvalue), val(signal_version)

    output:
    path "signalp_out"
    val pvalue
    val signal_version

    script:
    """
    signalp6 --organism ${organism} --fastafile ${fasta} --output_dir signalp_out --mode ${mode}
    """
}


process SIGNALP_PARSER {
    input:
    path out
    val pvalue
    val signalp_version
    val tsv_pro

    output:
    path "signalp_parsed.json"

    script:
    """
    python3 $projectDir/scripts/signalp/parser.py \
        ${out}/prediction_results.txt \
        ${pvalue} \
        ${signalp_version} \
        > signalp_parsed.json
    """
}
