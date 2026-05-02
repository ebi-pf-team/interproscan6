include { PREFILTER_SMART; PREPARE_SMART; SEARCH_SMART; PARSE_SMART } from  "../../modules/smart"

workflow SMART {
    take:
    ch_fasta     // [meta, fasta]
    ch_hmm      // [hmm3, hmm2]

    main:
    ch_hmm3 = ch_hmm.map { hmm3, hmm2 -> hmm3 }
    ch_hmm2 = ch_hmm.map { hmm3, hmm2 -> hmm2 }

    PREFILTER_SMART(
        ch_fasta,
        ch_hmm3
    )

    PREPARE_SMART(
        PREFILTER_SMART.out
    )

    SEARCH_SMART(
        PREPARE_SMART.out,
        ch_hmm2
    )

    PARSE_SMART(
        SEARCH_SMART.out
    )

    emit:
    PARSE_SMART.out  // [ meta, json ]
}