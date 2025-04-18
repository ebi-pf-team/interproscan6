include { RUN_HMMER as SEARCH_NCBIFAM } from "../../../modules/hmmer"
include { PARSE_NCBIFAM               } from "../../../modules/ncbifam"

workflow NCBIFAM {
    take:
    ch_seqs             // channel of tuples (index, fasta file)
    dir
    hmm

    main:
    SEARCH_NCBIFAM(
        ch_seqs,
        dir,
        hmm,
        "-Z 61295632 --cut_tc"
    )

    ch_ncbifam = PARSE_NCBIFAM(
        SEARCH_NCBIFAM.out
    )

    emit:
    ch_ncbifam
}