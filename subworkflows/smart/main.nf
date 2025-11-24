include { PREFILTER_SMART; PREPARE_SMART; SEARCH_SMART; PARSE_SMART } from  "../../modules/smart"

workflow SMART {
    take:
    ch_seqs     // channel of tuples (index, fasta file)
    dirpath     // data directory path
    hmmer3_hmm  // path to the HMM3 file in the data dir -> datadir/hmmer3_hmm
    hmmer2_dir  // path to the HMM2 directory in the data dir -> datadir/hmmer2_dir
    batch_size  // number of sequences per sub-batch for searching

    main:
    PREFILTER_SMART(
        ch_seqs,
        dirpath,
        hmmer3_hmm
    )

    PREPARE_SMART(
        PREFILTER_SMART.out,
        dirpath,
        hmmer2_dir,
        batch_size
    )

    SEARCH_SMART(
        PREPARE_SMART.out,
        dirpath,
        hmmer2_dir
    )

    ch_smart = PARSE_SMART(
        SEARCH_SMART.out,
        dirpath,
        hmmer2_dir
    )

    emit:
    ch_smart
}