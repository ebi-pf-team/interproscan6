include { RUN_PSSCAN; PARSE_PSSCAN } from  "../../../modules/pftools"

workflow PROSITE_PATTERNS {
    take:
    fasta       // [meta, fasta]
    files       // [dat, eval]

    main:
    RUN_PSSCAN(fasta, files)
    PARSE_PSSCAN(RUN_PSSCAN.out)

    emit:
    PARSE_PSSCAN.out  // [ meta, json ]
}