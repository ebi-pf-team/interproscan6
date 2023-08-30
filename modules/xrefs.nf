
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
    python $projectDir/scripts/xrefs.py -matches ${match_results} -entries ${entries} -go ${goterms} -pa ${pathways} > xrefs_results
    """
}