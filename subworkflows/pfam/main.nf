include { RUN_HMMER as SEARCH_PFAM } from "../../modules/hmmer"
include { PARSE_PFAM               } from  "../../modules/pfam"

workflow PFAM {
    take:
    ch_seqs          // channel of tuples (index, fasta file)
    dir              // str repr of the data directory path
    hmm              // str repr of the path to the HMM file in the data dir -> datadir/hmmFile
    dat              // str repr of the path to the DAT file in the data dir -> datadir/datFile

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