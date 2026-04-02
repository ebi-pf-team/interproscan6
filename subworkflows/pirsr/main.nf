include { RUN_HMMER as SEARCH_PIRSR } from  "../../modules/hmmer"
include { PARSE_PIRSR               } from  "../../modules/pirsr"

workflow PIRSR {
    take:
    fasta    // [meta, fasta]
    pirsr    // [hmm, json]

    main:
    hmm  = pirsr.map { hmm, json -> hmm }
    json = pirsr.map { hmm, json -> json }

    SEARCH_PIRSR(
        fasta,
        hmm,
        "-E 0.01 --acc"
    )

    PARSE_PIRSR(
        SEARCH_PIRSR.out,
        json
    )

    emit:
    PARSE_PIRSR.out  // [ meta, json ]
}