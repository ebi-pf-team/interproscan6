include { LOOKUP_CHECK } from "$projectDir/modules/local/lookup_check/main"
include { LOOKUP_MATCHES } from "$projectDir/modules/local/lookup_matches/main"
include { LOOKUP_NO_MATCHES } from "$projectDir/modules/local/lookup_no_matches/main"

workflow SEQUENCE_PRECALC {
    take:
    hash_sequence
    applications

    main:
    LOOKUP_CHECK(hash_sequence)
    LOOKUP_MATCHES(LOOKUP_CHECK.out, applications)
    LOOKUP_NO_MATCHES(LOOKUP_CHECK.out)

    emit:
    sequences_to_analyse = LOOKUP_NO_MATCHES.out
    parsed_matches = LOOKUP_MATCHES.out
}
