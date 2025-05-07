include { COMBINE_MATCHES               } from "../../modules/combine"
include { XREFS                         } from "../../modules/xrefs"
include { REPRESENTATIVE_LOCATIONS      } from "../../modules/representative_locations"

workflow COMBINE {
    take:
    match_results
    db_releases
    add_goterms
    add_pathways
    panther_paint_dir
    skip_interpro      // boolean used in production

    main:
    parsed_matches = Channel.empty()

    if (skip_interpro) {
        // COMBINE_MATCHES: Aggregate matches across all members for each sequence -> single JSON with all matches for the batch
        parsed_matches = COMBINE_MATCHES(match_results)
    } else {
        COMBINE_MATCHES(match_results)

        /* XREFS:
        Add signature and entry desc and names
        Add PAINT annotations (if panther is enabled)
        Add go terms (if enabled)
        Add pathways (if enabled)
        */
        XREFS(
            COMBINE_MATCHES.out,
            db_releases,
            add_goterms,
            add_pathways,
            panther_paint_dir
        )

        parsed_matches = REPRESENTATIVE_LOCATIONS(XREFS.out)
    }

    // Collect all JSON files into a single channel so we don't have cocurrent writing to the output files
    ch_results = parsed_matches
        .map { meta, json -> json }
        .collect()

    emit:
    ch_results
}