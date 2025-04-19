include { SEARCH_SFLD; POST_PROCESS_SFLD; PARSE_SFLD } from  "../../../modules/sfld"

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
        hmmfile
    )

    POST_PROCESS_SFLD(
        SEARCH_SFLD.out,
        dirpath,
        annofile
    )

    ch_sfld = PARSE_SFLD(
        POST_PROCESS_SFLD.out,
        dirpath,
        hierarchyfile
    )

    emit:
    ch_sfld
}