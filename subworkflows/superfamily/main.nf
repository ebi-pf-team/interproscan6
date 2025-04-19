include { SEARCH_SUPERFAMILY; PARSE_SUPERFAMILY } from  "../../modules/superfamily"

workflow SUPERFAMILY {
    take:
    ch_seqs
    dirpath
    hmm
    selfhits
    cla
    model
    pdbj95d

    main:
    SEARCH_SUPERFAMILY(
        ch_seqs,
        dirpath,
        hmm,
        selfhits,
        cla,
        model,
        pdbj95d
    )

    ch_superfams = PARSE_SUPERFAMILY(
        SEARCH_SUPERFAMILY.out,
        dirpath,
        model,
        hmm
    )

    emit:
    ch_superfams
}