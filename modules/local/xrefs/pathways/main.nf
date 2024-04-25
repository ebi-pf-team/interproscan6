process PATHWAYS {
    input:
    path matches
    val pathways

    output:
    path "xrefs_pathways_${matches}"

    script:
    """
    python3 $projectDir/scripts/xrefs/pathways.py ${matches} ${pathways} > xrefs_pathways_${matches}
    """
}
