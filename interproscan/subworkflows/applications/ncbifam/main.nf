include { RUN_HMMER     } from "../../../modules/hmmer"
include { PARSE_NCBIFAM } from "../../../modules/ncbifam"

workflow NCBIFAM {
    take:
    ch_seqs             // channel of tuples (index, fasta file)
    hmm                 // path to hmm file

    main:
    RUN_HMMER(
        ch_seqs,
        hmm,
        "-Z 61295632 --cut_tc"
    )

    ch_ncbifam = PARSE_NCBIFAM(
        RUN_HMMER.out
    )

    emit:
    ch_ncbifam
}