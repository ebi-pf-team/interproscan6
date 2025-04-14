include { SEARCH_SFLD; POST_PROCESS_SFLD; PARSE_SFLD } from  "../../../modules/sfld"

workflow SFLD {
    take:
    ch_seqs
    sfld_hmm
    sites_annotation
    sfld_hierarchy

    main:
    SEARCH_SFLD(
        ch_seqs,
        sfld_hmm
    )

    POST_PROCESS_SFLD(
        SEARCH_SFLD.out,
        sites_annotation
    )

    ch_sfld = PARSE_SFLD(
        POST_PROCESS_SFLD.out,
        sfld_hierarchy
    )

    emit:
    ch_sfld
}