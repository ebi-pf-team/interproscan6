process PATHWAYS {
    input:
    path matches
    val pathways

    output:
    path "pathways_${matches}"

    script:
    """
    python3 $projectDir/scripts/xrefs/pathways.py ${matches} ${pathways} > pathways_${matches}
    """
}
