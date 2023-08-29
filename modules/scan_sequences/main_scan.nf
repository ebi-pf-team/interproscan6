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

    sequences_application.map { sequences, application ->
        tuple(sequences, application)
    }
    .set{params_parser}


    HMMER_RUNNER(params_hmmer)
    PARSER(params_parser, HMMER_RUNNER.out)

    emit:
      PARSER.out
}
