include { HMMER_RUNNER } from "$projectDir/modules/scan_sequences/hmmer_runner"
include { PARSER } from "$projectDir/modules/scan_sequences/parser"

workflow MAIN_SCAN {
    take:
    fasta_application

    main:
    fasta_application.map { fasta, appl ->
        tuple(fasta, "${params.members_hmm[appl]}", "${params.members_switches[appl]}")
    }
    .set{hmmer_params}

    HMMER_RUNNER(hmmer_params)
    PARSER(HMMER_RUNNER.out)

    emit:
      PARSER.out
}
