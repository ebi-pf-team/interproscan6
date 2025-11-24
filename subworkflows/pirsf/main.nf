include { SEARCH_PIRSF; PARSE_PIRSF } from  "../../modules/pirsf"

workflow PIRSF {
    take:
    ch_seqs           // channel of tuples (index, fasta file)
    dirpath           // str repr of the data directory path
    hmmfile           // str repr of the path to the HMM file in the data dir -> datadir/hmmFile
    datfile           // str repr of the path to the DAT file in the data dir -> datadir/datFile

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