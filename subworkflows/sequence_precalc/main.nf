include { LOOKUP_CHECK } from "$projectDir/modules/local/lookup_check/main"
include { LOOKUP_MATCHES } from "$projectDir/modules/local/lookup_matches/main"
include { REVERSE_PARSE_SEQUENCE } from "$projectDir/modules/local/reverse_parse_sequence/main"

workflow SEQUENCE_PRECALC {
    take:
    hash_sequence
    applications

    main:
    LOOKUP_CHECK(hash_sequence)
    LOOKUP_MATCHES(LOOKUP_CHECK.out, applications)
    REVERSE_PARSE_SEQUENCE(LOOKUP_CHECK.out, hash_sequence)

    emit:
    sequences_to_analyse = REVERSE_PARSE_SEQUENCE.out
    parsed_matches = LOOKUP_MATCHES.out
}
