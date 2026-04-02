include { RUN_HMMER as SEARCH_ANTIFAM } from "../../modules/hmmer"
include { PARSE_ANTIFAM               } from "../../modules/antifam"

workflow ANTIFAM {
    take:
    ch_seqs             // channel of tuples (index, fasta file)
    hmmfile             // path to HMM

    main:
    SEARCH_ANTIFAM(
        ch_seqs,
        hmmfile,
        "--cut_ga"
    )

    PARSE_ANTIFAM(
        SEARCH_ANTIFAM.out
    )

    emit:
    PARSE_ANTIFAM.out  // [ meta, json ]
}