process SIGNALP_PARSER {
    input:
    path out
    val tsv_pro

    output:
    path "signalp_parsed.json"

    script:
    """
    python3 $projectDir/scripts/signalp/parser.py ${out}/prediction_results.txt > signalp_parsed.json
    """
}
