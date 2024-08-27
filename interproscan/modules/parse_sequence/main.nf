process PARSE_SEQUENCE {
    label 'io'
    errorStrategy { task.exitStatus in [143,1] ? 'terminate' : 'ignore' }

    input:
    val fasta_file
    val original_fasta_file
    val nucleic
    /*
    When nucleic is true, $fasta_file will contain the translated protein
    sequences from the predicted ORFs. So $originl_fasta_file
    is needed so that original nucleotide sequences can be
    associated with the corresponding ORF in the
    final output.
    */
    val applications

    output:
    path "parsed_sequences"

    script:
    """
    python3 $projectDir/interproscan/scripts/parse_sequence/parse_sequence.py \\
        ${fasta_file} \\
        ${original_fasta_file} \\
        ${nucleic} \\
        ${applications} \\
        "parsed_sequences"
    """
}
