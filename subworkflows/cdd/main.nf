include { RUN_RPSBLAST; RUN_RPSPROC; PARSE_RPSPROC } from  "../../modules/cdd"

workflow CDD {
    take:
    ch_seqs
    cdd_dir
    rpsblast_db
    rpsproc_db

    main:
    RUN_RPSBLAST(
        ch_seqs,
        cdd_dir,
        rpsblast_db
    )

    RUN_RPSPROC(
        RUN_RPSBLAST.out,
        cdd_dir,
        rpsproc_db
    )

    ch_cdd = PARSE_RPSPROC(RUN_RPSPROC.out)

    emit:
    ch_cdd
}