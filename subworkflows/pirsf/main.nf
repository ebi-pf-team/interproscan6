include { SEARCH_PIRSF; PARSE_PIRSF } from  "../../modules/pirsf"

workflow PIRSF {
    take:
    ch_seqs
    dirpath
    hmmfile
    datfile

    main:
    SEARCH_PIRSF(
        ch_seqs,
        dirpath,
        hmmfile
    )

    ch_pirsf = PARSE_PIRSF(
        SEARCH_PIRSF.out,
        dirpath,
        datfile
    )

    emit:
    ch_pirsf
}