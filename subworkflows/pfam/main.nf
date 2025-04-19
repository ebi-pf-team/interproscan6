include { RUN_HMMER as SEARCH_PFAM } from "../../modules/hmmer"
include { PARSE_PFAM               } from  "../../modules/pfam"

workflow PFAM {
    take:
    ch_seqs
    dir
    hmm
    dat

    main:
    SEARCH_PFAM(
        ch_seqs,
        dir,
        hmm,
        "-Z 61295632 --cut_ga"
    )

    ch_pfam = PARSE_PFAM(
        SEARCH_PFAM.out,
        dir,
        dat
    )

    emit:
    ch_pfam
}