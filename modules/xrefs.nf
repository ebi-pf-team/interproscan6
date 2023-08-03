process XREFS {
    input:
    val match_lookup
    val entries
    val goterms
    val pathways

    output:
    path 'match_lookup_xrefs'

    script:
    """
    python $projectDir/scripts/lookup/xrefs.py -matches ${match_lookup} -entries ${entries} -go ${goterms} -pa ${pathways} > match_lookup_xrefs
    """
}