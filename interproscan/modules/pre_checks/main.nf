process CHECK_NUCLEIC {
    label 'io'

    input:
    path fasta_file

    output:
    val ""  // need so the linter do not complain

    script:
    """
    python3 $projectDir/scripts/pre_checks/check_nucleic_seq.py ${fasta_file}
    """
}