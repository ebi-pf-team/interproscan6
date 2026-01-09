include { RUN_HAMAP; PARSE_HAMAP } from  "../../modules/hamap"

workflow HAMAP {
    take:
    ch_seqs        // channel of tuples (index, fasta file)
    hamap_dir      // str repr of the data directory path
    profiles_dir   // str repr of the path to the profiles directory in the data dir -> datadir/profiles

    main:
    RUN_HAMAP(
        ch_seqs,
        hamap_dir,
        profiles_dir
    )

    ch_hamap = PARSE_HAMAP(RUN_HAMAP.out)

    emit:
    ch_hamap
}