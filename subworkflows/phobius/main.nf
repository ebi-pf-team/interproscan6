include { WRITE_FASTA; SEARCH_PHOBIUS; PARSE_PHOBIUS } from  "../../modules/phobius"

workflow PHOBIUS {
    take:
    ch_seqs      // channel of tuples (index, fasta file)
    phobius_dir  // str repr of the data directory path

    main:
    WRITE_FASTA(ch_seqs)

    SEARCH_PHOBIUS(
        WRITE_FASTA.out,
        phobius_dir
    )

    ch_phobius = PARSE_PHOBIUS(SEARCH_PHOBIUS.out)

    emit:
    ch_phobius
}