process LOOKUP_NO_MATCHES {
    input:
    val checked_lookup

    output:
    path "no_match_lookup_fasta.fasta", optional: true

    script:
    """
    output=\$(python3 $projectDir/scripts/lookup/lookup_no_matches.py -checked "${checked_lookup}")

    if [ -n "\$output" ]; then
        echo "\$output" > no_match_lookup_fasta.fasta
    fi
    """
}
