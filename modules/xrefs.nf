process XREFS {
    input:
    val hash_seq
    val matches
    val entries
    val goterms
    val pathways

    output:
    path 'xrefs_results'

    script:
    """
    python3 $projectDir/scripts/xrefs.py -seq ${hash_seq} -matches ${matches} -entries ${entries} ${goterms ? "-go ${goterms}" : ""} ${pathways ? "-pa ${pathways}" : ""} > xrefs_results
    """
}