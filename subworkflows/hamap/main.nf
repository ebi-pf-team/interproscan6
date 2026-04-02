include { RUN_PFSEARCH ; PARSE_PFSEARCH } from  "../../modules/pftools"

workflow HAMAP {
    take:
    fasta          // [meta, fasta]
    profiles_dir   // path to the profiles directory

    main:
    RUN_PFSEARCH(
        fasta,
        profiles_dir
    )

    ch_hamap = PARSE_PFSEARCH(
        RUN_PFSEARCH.out,
        "HAMAP",
        "",
    )

    emit:
    ch_hamap
}
