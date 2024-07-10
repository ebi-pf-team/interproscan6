process PATHWAYS {
    label 'io'
    
    input:
    path matches
    val pathways

    output:
    path "pathways_${matches}"

    script:
    """
    python3 $projectDir/interproscan/scripts/xrefs/pathways.py ${matches} ${pathways} > pathways_${matches}
    """
}
