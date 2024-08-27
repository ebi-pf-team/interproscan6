process GOTERMS {
    label 'io'

    input:
    path matches
    val goterms

    output:
    path "goterms_${matches}"

    script:
    """
    python3 $projectDir/interproscan/scripts/xrefs/goterms.py \\
        ${matches} \\
        ${goterms} \\
        goterms_${matches}
    """
}
