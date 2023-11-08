include { LOOKUP_CHECK } from "$projectDir/modules/local/lookup_check/main"
include { LOOKUP_MATCHES } from "$projectDir/modules/local/lookup_matches/main"

workflow SEQUENCE_PRECALC {
    take:
    hash_sequence
    applications

    main:
    LOOKUP_CHECK(hash_sequence)
    LOOKUP_MATCHES(LOOKUP_CHECK.out, applications)

    emit:
    sequences_to_analyse = LOOKUP_CHECK.out
    parsed_matches = LOOKUP_MATCHES.out
}
