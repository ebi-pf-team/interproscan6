include { HMMER_RUNNER } from "$projectDir/modules/local/hmmer/runner/main"
include { HMMER_PARSER } from "$projectDir/modules/local/hmmer/parser/main"
include { GENE3D_POST_PROCESSER; FUNFAM_POST_PROCESSER; SFLD_POST_PROCESSER } from "$projectDir/modules/local/hmmer/post_processing/main"

workflow SEQUENCE_ANALYSIS {
    take:
    fasta
    applications

    main:
    Channel.from(applications.split(','))
    .branch { member ->
        def runner = ''
        if (params.members."${member}".runner == "hmmer") {
            runner = 'hmmer'
        }

        hmmer: runner == 'hmmer'
            return [ params.members."${member}".data, params.members."${member}".switches, member in ['funfam', 'gene3d', 'sfld'] ]
            /*
            The post processing of SFLD, FunFam and Gene3D HMMER hits requires the alignment file
            But only generate alignmnets for these tool to reduce volume and computational requirements
            */
        other: true
            log.info "Application ${member} (still) not supported"
    }.set { member_params }

    runner_hmmer_params = fasta.combine(member_params.hmmer)
    HMMER_RUNNER(runner_hmmer_params)
    HMMER_PARSER(HMMER_RUNNER.out, params.tsv_pro)

    emit:
    HMMER_PARSER.out
}
