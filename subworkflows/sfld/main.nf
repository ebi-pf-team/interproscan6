include { SEARCH_SFLD; PARSE_SFLD } from  "../../modules/sfld"

workflow SFLD {
    take:
    ch_seqs       // channel of tuples (index, fasta file)
    dirpath       // str repr of the data directory path
    hmmfile       // str repr of the path to the HMM file in the data dir -> datadir/hmmFile
    annofile      // str repr of the path to the annotation file in the data dir -> datadir/annoFile
    hierarchyfile // str repr of the path to the hierarchy file in the data dir -> datadir/hierarchyFile

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