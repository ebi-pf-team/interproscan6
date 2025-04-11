include { SEARCH_PIRSF; PARSE_PIRSF } from  "../../../modules/pirsf"

workflow PIRSF {
    take:
    ch_seqs
    pirsf_hmm
    pirsf_dat

    main:
    SEARCH_PIRSF(
        ch_seqs,
        pirsf_hmm
    )

    ch_pirsf = PARSE_PIRSF(
        SEARCH_PIRSF.out,
        pirsf_dat
    )

    emit:
    ch_pirsf
}