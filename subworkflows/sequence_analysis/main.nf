include { HMMER_RUNNER } from "$projectDir/modules/local/hmmer_runner/main"
include { HMMER_PARSER } from "$projectDir/modules/local/hmmer_parser/main"

workflow SEQUENCE_ANALYSIS {
    take:
    fasta_application

    main:
        fasta_application.map { fasta, appl ->
            runner: $params.members[appl].runner
            fasta_path: fasta
            member_configs: $params.members[appl]
        }.set {runner_params}

       if (runner_params.runner == "hmmer") {
            HMMER_RUNNER(runner_params.fasta_path, runner_params.member_configs)
            HMMER_PARSER(HMMER_RUNNER.out)
            result = HMMER_PARSER.out
        }
        else {
            error "Runner not supported"
        }

    emit:
      result
}
