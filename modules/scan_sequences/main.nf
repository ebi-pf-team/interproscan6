include { HMMER_RUNNER } from "$projectDir/modules/scan_sequences/hmmer_runner"
include { PARSER } from "$projectDir/modules/scan_sequences/parser"

workflow MAIN_SCAN {
    take:
    sequences_application

    main:
    sequences_application.map { sequences, application ->
        tuple(sequences, "${params.members_hmm[application]}", "${params.members_switches[application]}")
    }
    .set{params_hmmer}

    HMMER_RUNNER(params_hmmer)
//     PARSER(sequences_appl, HMMER_RUNNER.out)
}
