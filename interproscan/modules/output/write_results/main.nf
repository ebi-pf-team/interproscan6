process WRITE_RESULTS {
    label 'write_results'

    input:
    val file_name
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
    python3 -c "import json, sys; d = {}; [d.update(json.load(open(f))) for f in sys.argv[1:]]; print(json.dumps(d))" \\
     ${sequences.join(" ")} > sequences_hash.json

    python3 $projectDir/interproscan/scripts/output/write_results.py \\
        sequences_hash.json \\
        ${matches} \\
        ${format} \\
        ${file_name} \\
        $version \\
        ${nucleic} > debug_out
    """
}
