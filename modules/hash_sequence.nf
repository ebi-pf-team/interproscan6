process XREFS {
    input:
    val match_results
    val entries
    val goterms
    val pathways

    output:
    path 'xrefs_results'

    script:
    """
    python3 $projectDir/scripts/xrefs.py -matches ${match_results} -entries ${entries} ${goterms ? "-go ${goterms}" : ""} ${pathways ? "-pa ${pathways}" : ""} > xrefs_results
    """
}