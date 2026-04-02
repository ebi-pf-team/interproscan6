include { RUN_PFSEARCH ; PARSE_PFSEARCH } from  "../../../modules/pftools"

workflow PROSITE_PROFILES {
    take:
    fasta       // [meta, fasta]
    profiles    // [dir, skipped]

    main:
    ch_dir = profiles. map { dir, skipped -> dir }
    ch_skipped = profiles. map { dir, skipped -> skipped }

    RUN_PFSEARCH(
        fasta,
        ch_dir
    )

    PARSE_PFSEARCH(
        RUN_PFSEARCH.out,
        "PROSITE profiles",
        ch_skipped
    )

    emit:
    PARSE_PFSEARCH.out
}
