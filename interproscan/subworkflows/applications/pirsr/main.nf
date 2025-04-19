include { RUN_HMMER as SEARCH_PIRSR } from  "../../../modules/hmmer"
include { PARSE_PIRSR               } from  "../../../modules/pirsr"

workflow PIRSR {
    take:
    ch_seqs
    dirpath
    hmmfile
    rulesfile

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