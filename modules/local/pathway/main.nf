process XREFS {
    input:
    path matches
    val pathways

    output:
    path "pathway_xrefs_results_${matches}"

    script:
    """
    python3 $projectDir/scripts/xrefs.py -matches ${matches} ${entries} > xrefs_results_${matches}
    """
}
