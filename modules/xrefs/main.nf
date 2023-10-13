process XREFS {
    input:
    path matches
    val entries
    val goterms
    val pathways

    output:
    path "xrefs_results_${matches}"

    script:
    """
    python3 $projectDir/scripts/xrefs.py -matches ${matches} -entries ${entries} ${goterms ? "-go ${goterms}" : ""} ${pathways ? "-pa ${pathways}" : ""} > xrefs_results_${matches}
    """
}