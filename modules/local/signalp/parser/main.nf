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
