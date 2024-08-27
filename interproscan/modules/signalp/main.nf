process SIGNALP_RUNNER {
    container 'docker.io/library/signalp6'
    label 'signalp_runner'

    input:
    tuple path(fasta), val(mode), val(organism), val(pvalue), val(signal_version)

    output:
    path "signalp_out"
    val pvalue
    val signal_version
    val organism

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
    val organism

    output:
    path "signalp_parsed.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/signalp/parser.py \\
        ${out}/output.gff3 \\
        ${out}/prediction_results.txt \\
        ${pvalue} \\
        ${signalp_version} \\
        ${organism} \\
        signalp_parsed.json
    """
}
