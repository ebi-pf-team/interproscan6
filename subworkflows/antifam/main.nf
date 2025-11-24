include { RUN_HMMER as SEARCH_ANTIFAM } from "../../modules/hmmer"
include { PARSE_ANTIFAM               } from "../../modules/antifam"

workflow ANTIFAM {
    take:
    ch_seqs             // channel of tuples (index, fasta file)
    dirpath             // str repr of the data directory path
    hmmfile             // str repr of the path to the HMM file in the data dir -> datadir/hmmFile

    main:
    SEARCH_ANTIFAM(
        ch_seqs,
        dirpath,
        hmmfile,
        "--cut_ga"
    )

    ch_antifam = PARSE_ANTIFAM(
        SEARCH_ANTIFAM.out
    )

    emit:
    ch_antifam
}