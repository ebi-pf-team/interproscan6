include { RUN_HMMER as SEARCH_PIRSR } from  "../../modules/hmmer"
include { PARSE_PIRSR               } from  "../../modules/pirsr"

workflow PIRSR {
    take:
    ch_fasta    // [meta, fasta]
    ch_pirsr    // [hmm, json]

    main:
    ch_hmm  = ch_pirsr.map { hmm, json -> hmm }
    ch_json = ch_pirsr.map { hmm, json -> json }

    SEARCH_PIRSR(
        ch_fasta,
        ch_hmm,
        "-E 0.01 --acc"
    )

    PARSE_PIRSR(
        SEARCH_PIRSR.out,
        ch_json
    )

    emit:
    PARSE_PIRSR.out  // [ meta, json ]
}