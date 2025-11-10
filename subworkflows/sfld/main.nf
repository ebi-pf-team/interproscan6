include { SEARCH_SFLD; PARSE_SFLD } from  "../../modules/sfld"

workflow SFLD {
    take:
    ch_seqs
    dirpath
    hmmfile
    annofile
    hierarchyfile

    main:
    SEARCH_SFLD(
        ch_seqs,
        dirpath,
        hmmfile,
        annofile
    )

    ch_sfld = PARSE_SFLD(
        SEARCH_SFLD.out,
        dirpath,
        hierarchyfile
    )

    emit:
    ch_sfld
}