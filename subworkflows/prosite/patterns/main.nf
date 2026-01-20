include { RUN_PSSCAN; PARSE_PSSCAN } from  "../../../modules/pftools"

workflow PROSITE_PATTERNS {
    take:
    ch_seqs         // channel of tuples (index, fasta file)
    dirpath         // str repr of the data directory path
    datfile         // str repr of the path to the dat file in the data dir -> datadir/datfile
    evafile         // str repr of the path to the eva file in the data dir -> datadir/evafile

    main:
    RUN_PSSCAN(
        ch_seqs,
        dirpath,
        datfile,
        evafile
    )

    ch_prosite = PARSE_PSSCAN(RUN_PSSCAN.out)

    emit:
    ch_prosite
}