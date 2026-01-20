include { RUN_PFSEARCH ; PARSE_PFSEARCH } from  "../../modules/pftools"

workflow HAMAP {
    take:
    ch_seqs        // channel of tuples (index, fasta file)
    hamap_dir      // str repr of the data directory path
    profiles_dir   // str repr of the path to the profiles directory in the data dir -> datadir/profiles

    main:
    RUN_PFSEARCH(
        ch_seqs,
        hamap_dir,
        profiles_dir
    )

    ch_hamap = PARSE_PFSEARCH(
        RUN_PFSEARCH.out,
        "HAMAP",
        "",
        ""
    )

    emit:
    ch_hamap
}
