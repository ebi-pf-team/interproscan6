process COILS_RUNNER {
    label 'coils_runner'

    input:
    tuple path(fasta), val(release), val(switches)

    output:
    path "coil_out"
    val release

    script:
    """
    /opt/coils/ncoils ${switches} < ${fasta} > coil_out
    """
}


process COILS_PARSER {
    label 'analysis_parser'

    input:
    path out
    val release

    output:
    path "coils_parsed.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/coils/parser.py \
        ${out} \\
        ${release} \\
        coils_parsed.json
    """
}
