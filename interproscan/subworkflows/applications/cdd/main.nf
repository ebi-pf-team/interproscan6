include { RUN_RPSBLAST; RUN_RPSPROC; PARSE_RPSPROC } from  "../../../modules/cdd"

workflow CDD {
    take:
    ch_seqs
    rpsblast_db
    rpsproc_db

    main:
    RUN_RPSBLAST(
        ch_seqs,
        rpsblast_db
    )

    RUN_RPSPROC(
        RUN_RPSBLAST.out,
        rpsproc_db
    )

    ch_cdd = PARSE_RPSPROC(RUN_RPSPROC.out)

    emit:
    ch_cdd
}