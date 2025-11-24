include { RUN_HMMER as SEARCH_PIRSR } from  "../../modules/hmmer"
include { PARSE_PIRSR               } from  "../../modules/pirsr"

workflow PIRSR {
    take:
    ch_seqs           // channel of tuples (index, fasta file)
    dirpath           // data directory path
    hmmfile           // path to the HMM file in the data dir -> datadir/hmmFile
    rulesfile         // path to the rules file in the data dir -> datadir/rulesFile

    main:
    SEARCH_PIRSR(
        ch_seqs,
        dirpath,
        hmmfile,
        "-E 0.01 --acc"
    )

    ch_pirsr = PARSE_PIRSR(
        SEARCH_PIRSR.out,
        dirpath,
        rulesfile
    )

    emit:
    ch_pirsr
}