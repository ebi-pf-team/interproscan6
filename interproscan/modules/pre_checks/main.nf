process CHECK_SEQUENCES {
    label 'io'

    input:
    path fasta_file
    val applications
    val is_nucleic
    path outdir

    output:
    val ""  // needed so the linter do not complain

    script:
    """
    python3 $projectDir/interproscan/scripts/pre_checks/check_seqs.py \\
        ${fasta_file} \\
        ${applications} \\
        ${is_nucleic} \\
        ${outdir}
    """
}
