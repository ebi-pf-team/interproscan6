process ENTRIES {
    input:
    path matches
    val entries

    output:
    path "xrefs_entries_${matches}"

    script:
    """
    python3 $projectDir/scripts/xrefs/entries.py ${matches} ${entries} > xrefs_entries_${matches}
    """
}
