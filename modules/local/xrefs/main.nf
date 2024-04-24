process XREFS {
    input:
    path matches
    val entries

    output:
    path "xrefs_results_${matches}"

    script:
    """
    python3 $projectDir/scripts/xrefs.py ${matches} ${entries} > xrefs_results_${matches}
    """
}
