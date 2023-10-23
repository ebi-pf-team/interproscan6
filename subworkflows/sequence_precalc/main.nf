include { LOOKUP_CHECK } from "$projectDir/modules/lookup_check/main"
include { LOOKUP_MATCHES } from "$projectDir/modules/lookup_matches/main"

workflow SEQUENCE_PRECALC {
    take:
    fasta_application

    main:
    LOOKUP_CHECK(fasta_application.map { it.first() })
    LOOKUP_MATCHES(fasta_application.map { it.last() }, LOOKUP_CHECK.out.map { it.first() })

    emit:
      LOOKUP_CHECK.out.map { it.last() }
      LOOKUP_MATCHES.out
}
