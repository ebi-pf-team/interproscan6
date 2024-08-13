process PARSE_SEQUENCE {
    label 'io'
    errorStrategy { task.exitStatus in [143,1] ? 'terminate' : 'ignore' }

    input:
    val fasta_file
    val applications

    output:
    path "parsed_sequences"

    script:
    """
    python3 $projectDir/interproscan/scripts/parse_sequence.py ${fasta_file} ${applications} > parsed_sequences
    """
}