process GOTERMS {
    input:
    path matches
    val goterms

    output:
    path "xrefs_goterm_${matches}"

    script:
    """
    python3 $projectDir/scripts/xrefs/goterms.py ${matches} ${goterms} > xrefs_goterms_${matches}
    """
}
