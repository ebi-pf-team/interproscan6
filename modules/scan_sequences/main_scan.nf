include { HMMER_RUNNER } from "$projectDir/modules/scan_sequences/hmmer_runner"
include { PARSER } from "$projectDir/modules/scan_sequences/parser"

process MAIN_SCAN {
    input:
    path fasta_file
    val application

    output:
    path output_main_scan

    script:
    """
    echo ${fasta_file} ${application} > output_main_scan
    """
}

workflow {
    hmm_app = params.members_hmm.hmm_${application}_path
    options = params.switches.switches_${application}

    HMMER_RUNNER(fasta_file, hmm_app, options)
    PARSER(HMMER_RUNNER.out, fasta_file)
}
