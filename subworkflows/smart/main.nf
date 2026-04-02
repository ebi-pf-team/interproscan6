include { PREFILTER_SMART; PREPARE_SMART; SEARCH_SMART; PARSE_SMART } from  "../../modules/smart"

workflow SMART {
    take:
    fasta       // [meta, fasta]
    hmm         // [hmm3, hmm2]

    main:
    hmm3 = hmm.map { hmm3, hmm2 -> hmm3 }
    hmm2 = hmm.map { hmm3, hmm2 -> hmm2 }

    PREFILTER_SMART(
        fasta,
        hmm3
    )

    PREPARE_SMART(
        PREFILTER_SMART.out,
        hmm2,
    )

    SEARCH_SMART(
        PREPARE_SMART.out,
        hmm2
    )

    PARSE_SMART(
        SEARCH_SMART.out,
        hmm2
    )

    emit:
    PARSE_SMART.out  // [ meta, json ]
}