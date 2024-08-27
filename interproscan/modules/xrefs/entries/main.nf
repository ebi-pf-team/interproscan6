process ENTRIES {
    label 'io'

    input:
    path matches
    val entries

    output:
    path "entries_${matches}"

    script:
    """
    python3 $projectDir/interproscan/scripts/xrefs/entries.py \\
        ${matches} \\
        ${entries} \\
        entries_${matches}
    """
}
