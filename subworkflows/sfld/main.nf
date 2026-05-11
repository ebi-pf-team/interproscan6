include { SEARCH_SFLD; PARSE_SFLD } from  "../../modules/sfld"

workflow SFLD {
    take:
    fasta       // [meta, fasta]
    files       // [hmm, sites, hierarchy]

    main:
    ch_search = files.map { hmm, sites, hierarchy -> tuple(hmm, sites) }
    ch_parse = files.map { hmm, sites, hierarchy -> hierarchy }

    SEARCH_SFLD(fasta, ch_search)

    PARSE_SFLD(SEARCH_SFLD.out, ch_parse)

    emit:
    PARSE_SFLD.out  // [ meta, json ]
}