include { PREFILTER_SMART; PREPARE_SMART; SEARCH_SMART; PARSE_SMART } from  "../../../modules/smart"

workflow SMART {
    take:
    ch_seqs
    hmmer3_hmm
    hmmer2_hmm
    smart_hmm_dir
    smart_chunksize

    main:
    PREFILTER_SMART(
        ch_seqs,
        hmmer3_hmm
    )

    PREPARE_SMART(
        PREFILTER_SMART.out,
        smart_chunksize,
        smart_hmm_dir
    )

    SEARCH_SMART(
        PREPARE_SMART.out,
        smart_hmm_dir
    )

    ch_smart = PARSE_SMART(
        SEARCH_SMART.out,
        hmmer2_hmm
    )

    emit:
    ch_smart
}