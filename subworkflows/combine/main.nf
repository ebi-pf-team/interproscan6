include { COMBINE_MATCHES; COMBINE_MATCHES_LOCAL                   } from "${moduleDir}/../../modules/combine"
include { XREFS; XREFS_LOCAL                                       } from "${moduleDir}/../../modules/xrefs"
include { REPRESENTATIVE_LOCATIONS; REPRESENTATIVE_LOCATIONS_LOCAL } from "${moduleDir}/../../modules/representative_locations"

workflow COMBINE {
    take:
    match_results
    db_releases
    add_goterms
    add_pathways
    panther_paint_dir
    skip_interpro      // boolean used in production
    batch_size

    main:
    parsed_matches = Channel.empty()

    if (skip_interpro) {
        // COMBINE_MATCHES: Aggregate matches across all members for each sequence -> single JSON with all matches for the batch
        if (batch_size < 50000) {
            parsed_matches = COMBINE_MATCHES_LOCAL(match_results)
        } else {
            parsed_matches = COMBINE_MATCHES(match_results)
        }
    } else {
        if (batch_size < 50000) {
            COMBINE_MATCHES_LOCAL(match_results)

            XREFS_LOCAL(
                COMBINE_MATCHES_LOCAL.out,
                db_releases,
                add_goterms,
                add_pathways,
                panther_paint_dir
            )

            parsed_matches = REPRESENTATIVE_LOCATIONS_LOCAL(XREFS_LOCAL.out)
        } else {
            COMBINE_MATCHES(match_results)

            XREFS(
                COMBINE_MATCHES.out,
                db_releases,
                add_goterms,
                add_pathways,
                panther_paint_dir
            )

            parsed_matches = REPRESENTATIVE_LOCATIONS(XREFS.out)
        }
    }

    // Collect all JSON files into a single channel so we don't have cocurrent writing to the output files
    ch_results = parsed_matches
        .map { meta, json -> json }
        .collect()

    emit:
    ch_results
}