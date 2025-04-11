include { RUN_HMMER     } from "../../../modules/hmmer"
include { PARSE_ANTIFAM } from "../../../modules/antifam"

workflow ANTIFAM {
    take:
    ch_seqs             // channel of tuples (index, fasta file)
    hmm                 // path to hmm file

    main:
    RUN_HMMER(
        ch_seqs,
        hmm,
        "--cut_ga"
    )

    ch_antifam = PARSE_ANTIFAM(
        RUN_HMMER.out
    )

    emit:
    ch_antifam
}