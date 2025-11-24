include { RUN_PFSEARCH ; PARSE_PFSEARCH } from  "../../../modules/prosite/profiles"

workflow PROSITE_PROFILES {
    take:
    ch_seqs        // channel of tuples (index, fasta file)
    dirpath        // data directory path
    profiles_dir   // path to the profiles directory in the data dir -> datadir/profilesDir
    blacklist_file // path to the blacklist file in the data dir -> datadir/blacklistFile

    main:
    RUN_PFSEARCH(
        ch_seqs,
        dirpath,
        profiles_dir
    )
    ch_prosite = PARSE_PFSEARCH(
        RUN_PFSEARCH.out,
        dirpath,
        blacklist_file
    )

    emit:
    ch_prosite
}