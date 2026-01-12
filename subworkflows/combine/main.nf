include { COMBINE_MATCHES; COMBINE_MATCHES_LOCAL                   } from "../../modules/combine"
include { XREFS; XREFS_LOCAL                                       } from "../../modules/xrefs"
include { REPRESENTATIVE_LOCATIONS; REPRESENTATIVE_LOCATIONS_LOCAL } from "../../modules/representative_locations"

workflow COMBINE {
    take:
    match_results
    db_releases
    add_goterms
    add_pathways
    panther_paint_dir
    skip_repr_locations      // boolean used in production
    batch_size

    main:
    ch_xrefs = Channel.empty()
    ch_combined = Channel.empty()

    if (batch_size < 50000) {
        COMBINE_MATCHES_LOCAL(match_results)

        ch_xrefs = XREFS_LOCAL(
            COMBINE_MATCHES_LOCAL.out,
            db_releases,
            add_goterms,
            add_pathways,
            panther_paint_dir
        )
    } else {
        COMBINE_MATCHES(match_results)

        ch_xrefs = XREFS(
            COMBINE_MATCHES.out,
            db_releases,
            add_goterms,
            add_pathways,
            panther_paint_dir
        )
    }

    if (skip_repr_locations) {
        ch_combined = ch_xrefs
    } else if (batch_size < 50000) {
        ch_combined = REPRESENTATIVE_LOCATIONS_LOCAL(ch_xrefs)
    } else {
        ch_combined = REPRESENTATIVE_LOCATIONS(ch_xrefs)
    }

    // Collect all JSON files into a single channel so we don't have cocurrent writing to the output files
    ch_results = ch_combined
        .map { meta, json -> json }
        .collect()

    emit:
    ch_results
}