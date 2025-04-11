include { RUN_HMMER as SEARCH_ANTIFAM } from "../../../modules/hmmer"
include { PARSE_ANTIFAM               } from "../../../modules/antifam"

workflow ANTIFAM {
    take:
    ch_seqs             // channel of tuples (index, fasta file)
    hmm                 // path to hmm file

    main:
    SEARCH_ANTIFAM(
        ch_seqs,
        hmm,
        "--cut_ga"
    )

    ch_antifam = PARSE_ANTIFAM(
        SEARCH_ANTIFAM.out
    )

    emit:
    ch_antifam
}