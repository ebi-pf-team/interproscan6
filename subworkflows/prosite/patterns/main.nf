include { RUN_PFSCAN; PARSE_PFSCAN } from  "../../../modules/prosite/patterns"

workflow PROSITE_PATTERNS {
    take:
    ch_seqs         // channel of tuples (index, fasta file)
    dirpath         // str repr of the data directory path
    datfile         // str repr of the path to the dat file in the data dir -> datadir/datfile
    evafile         // str repr of the path to the eva file in the data dir -> datadir/evafile

    main:
    RUN_PFSCAN(
        ch_seqs,
        dirpath,
        datfile,
        evafile
    )

    ch_prosite = PARSE_PFSCAN(RUN_PFSCAN.out)

    emit:
    ch_prosite
}