include { HMMER_RUNNER } from "$projectDir/modules/hmmer_runner/main"
include { HMMER_PARSER } from "$projectDir/modules/hmmer_parser/main"

workflow SEQUENCE_ANALYSIS {
    take:
    fasta_application
    tsv_pro

    main:
    fasta_application.map { fasta, appl ->
        tuple(fasta, "${params.members_hmm[appl]}", "${params.members_switches[appl]}")
    }
    .set{hmmer_params}

    HMMER_RUNNER(hmmer_params)
    PARSER(HMMER_RUNNER.out, tsv_pro)

    emit:
      PARSER.out
}
