include { RUN_HMMER as SEARCH_NCBIFAM } from "${moduleDir}/../../modules/hmmer"
include { PARSE_NCBIFAM               } from "${moduleDir}/../../modules/ncbifam"

workflow NCBIFAM {
    take:
    ch_seqs             // channel of tuples (index, fasta file)
    dir                 // str repr of the data directory path
    hmm                 // str repr of the path to the HMM file in the data dir -> datadir/hmmFile

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