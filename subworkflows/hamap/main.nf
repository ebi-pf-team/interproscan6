include { PREPROCESS_HAMAP; PREPARE_HAMAP; RUN_HAMAP; PARSE_HAMAP } from  "../../modules/hamap"

workflow HAMAP {
    take:
    ch_seqs        // channel of tuples (index, fasta file)
    hamap_dir      // str repr of the data directory path
    hmm_file       // str repr of the path to the HMM file in the data dir -> datadir/hmmFile
    profiles_dir   // str repr of the path to the profiles directory in the data dir -> datadir/profiles

    main:
    PREPROCESS_HAMAP(
        ch_seqs,
        hamap_dir,
        hmm_file
    )

    PREPARE_HAMAP(
        PREPROCESS_HAMAP.out,
        hamap_dir,
        profiles_dir
    )

    RUN_HAMAP(
        PREPARE_HAMAP.out,
        hamap_dir,
        profiles_dir
    )

    ch_hamap = PARSE_HAMAP(RUN_HAMAP.out)

    emit:
    ch_hamap
}