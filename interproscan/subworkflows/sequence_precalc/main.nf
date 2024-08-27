include {
    LOOKUP_CHECK;
    LOOKUP_MATCHES;
    LOOKUP_NO_MATCHES;
} from "$projectDir/interproscan/modules/lookup/main"

workflow SEQUENCE_PRECALC {
    take:
    hash_sequence
    applications
    is_test  // boolean, true = unit test so don't make calls to remote servers

    main:
    LOOKUP_CHECK(hash_sequence, is_test)
    LOOKUP_MATCHES(LOOKUP_CHECK.out, applications, is_test)
    LOOKUP_NO_MATCHES(LOOKUP_CHECK.out)

    if (LOOKUP_MATCHES.out.isEmpty()) {
        parsed_matches = null
    } else {
        parsed_matches = LOOKUP_MATCHES.out
    }

    emit:
    sequences_to_analyse = LOOKUP_NO_MATCHES.out
    parsed_matches
}
