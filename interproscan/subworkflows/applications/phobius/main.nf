include { SEARCH_PHOBIUS; PARSE_PHOBIUS } from  "../../../modules/phobius"

workflow PHOBIUS {
    take:
    ch_seqs
    phobius_dir

    main:
    SEARCH_PHOBIUS(
        ch_seqs,
        phobius_dir
    )

    ch_phobius = PARSE_PHOBIUS(SEARCH_PHOBIUS.out)

    emit:
    ch_phobius
}