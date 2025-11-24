include { SEARCH_CDD; PARSE_CDD } from  "../../modules/cdd"

workflow CDD {
    take:
    ch_seqs               // channel of tuples (index, fasta file)
    cdd_dir               // data directory path
    rpsblast_db           // path to the RPS-BLAST database in the data dir -> datadir/rpsblastDB
    rpsproc_db            // path to the RPS-PROC database in the data dir -> datadir/rpsprocDB

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