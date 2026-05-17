include { SEARCH_PIRSF; PARSE_PIRSF } from  "../../modules/pirsf"

workflow PIRSF {
    take:
    fasta    // [meta, fasta]
    pirsf    // [hmm, dat]

    main:
    ch_hmm = pirsf.map { hmm, _dat -> hmm }
    ch_dat = pirsf.map { _hmm, dat -> dat }

    SEARCH_PIRSF(
        fasta,
        ch_hmm
    )

    PARSE_PIRSF(
        SEARCH_PIRSF.out,
        ch_dat
    )

    emit:
    PARSE_PIRSF.out  // [ meta, json ]
}
