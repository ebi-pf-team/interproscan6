process COILS_RUNNER {
    container 'docker.io/library/coils'
    label 'coils_runner'

    input:
    tuple path(fasta), val(release) val(switches)

    output:
    path "coils_out.txt"
    val version

    script:
    """
    ncoils ${switches} < ${fasta} > coils_out.txt
    """
}


process COILS_PARSER {
    label 'analysis_parser'

    input:
    path out
    val version

    output:
    path "coils_parsed.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/coils/parser.py \
        ${out} \
        ${version} \
        > coils_parsed.json
    """
}
