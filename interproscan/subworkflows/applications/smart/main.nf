include { PREFILTER_SMART; PREPARE_SMART; SEARCH_SMART; PARSE_SMART } from  "../../../modules/smart"

workflow SMART {
    take:
    ch_seqs
    dirpath
    hmmer3_hmm
    hmmer2_dir
    chunksize

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
        chunksize
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