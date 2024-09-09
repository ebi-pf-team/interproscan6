process CHECK_NUCLEIC {
    label 'io'

    input:
    path fasta_file

    output:
    val ""  // needed so the linter do not complain

    script:
    """
    python3 $projectDir/interproscan/scripts/pre_checks/check_nucleic_seq.py ${fasta_file}
    """
}