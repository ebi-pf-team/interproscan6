process SIGNALP_PARSER {
    input:
    path out
    val tsv_pro

    output:
    path "signalp_parsed"

    script:
    """
    python3 $projectDir/scripts/signalp/???? > signalp_parsed.json
    """
}
