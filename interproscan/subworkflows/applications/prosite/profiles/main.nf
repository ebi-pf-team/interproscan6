include { RUN_PFSEARCH ; PARSE_PFSEARCH } from  "../../../../modules/prosite/profiles"

workflow PROSITE_PROFILES {
    take:
    ch_seqs
    prosite_profiles_dir
    skip_flagged_profiles

    main:
    RUN_PFSEARCH(
        ch_seqs,
        prosite_profiles_dir
    )
    ch_prosite = PARSE_PFSEARCH(
        RUN_PFSEARCH.out,
        skip_flagged_profiles
    )

    emit:
    ch_prosite
}