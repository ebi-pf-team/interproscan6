process GOTERMS {
    input:
    path matches
    val goterms

    output:
    path "goterm_${matches}"

    script:
    """
    python3 $projectDir/scripts/xrefs/goterms.py ${matches} ${goterms} > goterms_${matches}
    """
}
