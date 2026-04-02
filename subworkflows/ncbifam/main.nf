include { RUN_HMMER as SEARCH_NCBIFAM } from "../../modules/hmmer"
include { PARSE_NCBIFAM               } from "../../modules/ncbifam"

workflow NCBIFAM {
    take:
    fasta     // [meta, fasta]
    hmm       // path to the HMM file

    main:
    SEARCH_NCBIFAM(
        fasta,
        hmm,
        "-Z 61295632 --cut_tc"
    )

    ch_ncbifam = PARSE_NCBIFAM(
        SEARCH_NCBIFAM.out
    )

    emit:
    ch_ncbifam
}