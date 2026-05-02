include { RUN_HMMER as SEARCH_PFAM } from "../../modules/hmmer"
include { PARSE_PFAM               } from "../../modules/pfam"

workflow PFAM {
    take:
    ch_fasta  // [meta, fasta] 
    pfam     // [hmm, dat]

    main:
    ch_hmm = pfam.map { hmm, dat -> hmm }
    SEARCH_PFAM(
        ch_fasta,
        ch_hmm,
        "-Z 61295632 --cut_ga"
    )

    ch_dat = pfam.map { hmm, dat -> dat }
    PARSE_PFAM(
        SEARCH_PFAM.out,
        ch_dat
    )

    emit:
    PARSE_PFAM.out  // [ meta, json ]
}