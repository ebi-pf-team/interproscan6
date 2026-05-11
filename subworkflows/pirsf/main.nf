include { SEARCH_PIRSF; PARSE_PIRSF } from  "../../modules/pirsf"

workflow PIRSF {
    take:
    fasta    // [meta, fasta]
    pirsf    // [hmm, dat]

    main:
    hmm = pirsf.map { hmm, dat -> hmm }
    dat = pirsf.map { hmm, dat -> dat }

    SEARCH_PIRSF(
        fasta,
        hmm
    )

    PARSE_PIRSF(
        SEARCH_PIRSF.out,
        dat
    )

    emit:
    PARSE_PIRSF.out  // [ meta, json ]
}