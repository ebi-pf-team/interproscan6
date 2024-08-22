process WRITE_RESULTS {
    label 'write_results'

    input:
    val fasta
    val sequences
    path matches
    val format
    val output_path
    val version
    val nucleic

    output:
    path "*.ips6.*"

    script:
    """
    file_name=\$(basename $fasta)
    cat ${sequences.join(" ")} > sequences_hash.json
    python3 $projectDir/interproscan/scripts/output/write_results.py \\
        sequences_hash.json \\
        ${matches} \\
        ${format} \\
        \${file_name} \\
        $version \\
        ${nucleic} > debug_out
    """
}
