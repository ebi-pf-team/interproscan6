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

    PARSE_PFSEARCH(
        RUN_PFSEARCH.out,
        "HAMAP",
        "",
    )

    emit:
    PARSE_PFSEARCH.out  // [ meta, json ]
}
