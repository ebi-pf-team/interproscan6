include { SEARCH_CDD; PARSE_CDD } from  "../../modules/cdd"

workflow CDD {
    take:
    ch_seqs
    cdd_dir
    rpsblast_db
    rpsproc_db

    main:
    SEARCH_CDD(
        ch_seqs,
        cdd_dir,
        rpsblast_db,
        rpsproc_db
    )

    ch_cdd = PARSE_CDD(SEARCH_CDD.out)

    emit:
    ch_cdd
}