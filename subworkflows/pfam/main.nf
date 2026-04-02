include { RUN_HMMER as SEARCH_PFAM } from "../../modules/hmmer"
include { PARSE_PFAM               } from "../../modules/pfam"

workflow PFAM {
    take:
    fasta  // [meta, fasta] 
    pfam   // [hmm, dat]

    main:
    hmm = pfam.map { hmm, dat -> hmm }
    SEARCH_PFAM(
        fasta,
        hmm,
        "-Z 61295632 --cut_ga"
    )

    dat = pfam.map { hmm, dat -> dat }
    PARSE_PFAM(
        SEARCH_PFAM.out,
        dat
    )

    emit:
    PARSE_PFAM.out
}