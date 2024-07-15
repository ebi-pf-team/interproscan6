include { LOOKUP_CHECK } from "$projectDir/interproscan/modules/lookup/check/main"
include { LOOKUP_MATCHES } from "$projectDir/interproscan/modules/lookup/matches/main"
include { LOOKUP_NO_MATCHES } from "$projectDir/interproscan/modules/lookup/no_matches/main"

workflow SEQUENCE_PRECALC {
    take:
    hash_sequence
    applications
    is_test  // boolean, true = unit test so don't make calls to remote servers

    main:
    LOOKUP_CHECK(hash_sequence, is_test)
    LOOKUP_MATCHES(LOOKUP_CHECK.out, applications, is_test)
    LOOKUP_NO_MATCHES(LOOKUP_CHECK.out)

    emit:
    sequences_to_analyse = LOOKUP_NO_MATCHES.out
    parsed_matches = LOOKUP_MATCHES.out
}
