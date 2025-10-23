include { WRITE_FASTA; SEARCH_PHOBIUS; PARSE_PHOBIUS } from  "../../modules/phobius"

workflow PHOBIUS {
    take:
    ch_seqs
    phobius_dir

    main:
    ch_split = ch_seqs
        .splitFasta( by: 5000, file: true )

    WRITE_FASTA(ch_split)

    SEARCH_PHOBIUS(
        WRITE_FASTA.out,
        phobius_dir
    )

    ch_phobius = PARSE_PHOBIUS(SEARCH_PHOBIUS.out)

    emit:
    ch_phobius
}