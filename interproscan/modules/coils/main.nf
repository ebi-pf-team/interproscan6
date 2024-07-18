process COILS_RUNNER {
    container 'docker.io/biocontainers/ncoils:v2002-7-deb_cv1'
    label 'coils_runner'

    input:
    tuple path(fasta), val(release), val(switches)

    output:
    path "coil_out"
    val release

    script:
    """
    ncoils ${switches} < ${fasta} &> coil_out
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
