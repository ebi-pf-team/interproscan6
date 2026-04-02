include { SEARCH_CDD; PARSE_CDD } from  "../../modules/cdd"

workflow CDD {
    take:
    fasta       // [meta, fasta]
    rps         // [dir, blast_db, proc_db]

    main:
    SEARCH_CDD(
        fasta,
        rps
    )

    PARSE_CDD(SEARCH_CDD.out)

    emit:
    PARSE_CDD.out // [ meta, json ]
}