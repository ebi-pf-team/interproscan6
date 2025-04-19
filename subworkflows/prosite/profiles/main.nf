include { RUN_PFSEARCH ; PARSE_PFSEARCH } from  "../../../modules/prosite/profiles"

workflow PROSITE_PROFILES {
    take:
    ch_seqs
    dirpath
    profiles_dir
    blacklist_file

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