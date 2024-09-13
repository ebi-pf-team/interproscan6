process COILS_RUNNER {
    label 'coils_runner'

    input:
    tuple path(fasta), val(switches)

    output:
    path "coil_out"

    script:
    """
    /opt/coils/ncoils ${switches} < ${fasta} > coil_out
    """
}


process COILS_PARSER {
    label 'analysis_parser'

    input:
    path out

    output:
    path "coils_parsed.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/coils/parser.py \\
        ${out} \\
        coils_parsed.json
    """
}
