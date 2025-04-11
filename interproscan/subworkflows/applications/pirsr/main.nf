include { RUN_HMMER as SEARCH_PIRSR } from  "../../../modules/hmmer"
include { PARSE_PIRSR               } from  "../../../modules/pirsr"

workflow PIRSR {
    take:
    ch_seqs
    pirsr_hmm
    pirsr_rules

    main:
    SEARCH_PIRSR(
        ch_seqs,
        pirsr_hmm,
        "-E 0.01 --acc"
    )

    ch_pirsr = PARSE_PIRSR(
        SEARCH_PIRSR.out,
        pirsr_rules
    )

    emit:
    ch_pirsr
}