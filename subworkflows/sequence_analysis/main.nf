include { HMMER_RUNNER as GENERIC_HMMER_RUNNER; HMMER_RUNNER as SFLD_HMMER_RUNNER } from "$projectDir/modules/local/hmmer/runner/main"
include { HMMER_PARSER } from "$projectDir/modules/local/hmmer/parser/main"
include { GENE3D_POST_PROCESSER; FUNFAM_POST_PROCESSER; SFLD_POST_PROCESSER } from "$projectDir/modules/local/hmmer/post_processing/main"
include { SFLD_PARSER } from "$projectDir/modules/local/hmmer/parser/slfd"

workflow SEQUENCE_ANALYSIS {
    take:
    fasta
    applications

    main:
    // Divide members up into their respective analysis pipelines/methods
    Channel.from(applications.split(','))
    .branch { member ->
        def runner = ''
        if (params.members."${member}".runner == "hmmer") {
            runner = 'hmmer'
        }
        // funfam
        // gene3d
        // hamap
        // panther
        // pfam
        // pirsf
        // pirsr
        // smart ?
        // superfamily
        if (member == 'sfld') {
            runner = 'sfld'
        }

        /*
        The post processing of SFLD, FunFam and Gene3D HMMER hits requires additional files
        and parameters relative to the generic hmmer runner and parser in IPS6
        */
        hmmer: runner == 'hmmer'
            return [ 
                params.members."${member}".hmm, params.members."${member}".switches, 
                false, []
            ]
        sfld : runner == 'sfld'
            return [
                params.members."${member}".hmm, params.members."${member}".switches,
                true, [
                    params.members."${member}".postprocess.bin,
                    params.members."${member}".postprocess.sites_annotation,
                    params.members."${member}".postprocess.hierarchy
                ]
            ]
        other: true
            log.info "Application ${member} (still) not supported"
    }.set { member_params }

    runner_hmmer_params = fasta.combine(member_params.hmmer)
    GENERIC_HMMER_RUNNER(runner_hmmer_params)
    HMMER_PARSER(GENERIC_HMMER_RUNNER.out, params.tsv_pro)

    runner_hmmer_sfld_params = fasta.combine(member_params.sfld)
    SFLD_HMMER_RUNNER(runner_hmmer_sfld_params)
    SFLD_POST_PROCESSER(SFLD_HMMER_RUNNER.out, params.tsv_pro)
    // SFLD_PARSER(SFLD_POST_PROCESSER.out, member_params.sfld, params.tsv_pro)

    emit:
    "HMMER_PARSER.out"
}
