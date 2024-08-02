process SIGNALP_RUNNER {
    label 'signalp_runner'

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
    label 'analysis_parser'

    input:
    path out
    val pvalue
    val signalp_version

    output:
    path "signalp_parsed.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/signalp/parser.py \
        ${out}/prediction_results.txt \
        ${pvalue} \
        ${signalp_version} \
        > signalp_parsed.json
    """
}
