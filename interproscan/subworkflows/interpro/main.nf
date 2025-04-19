include { XREFS                         } from "../../modules/xrefs"
include { REPRESENTATIVE_LOCATIONS      } from "../../modules/representative_locations"

workflow INTERPRO {
    take:
    match_results
    applications
    db_releases
    add_goterms
    add_pathways
    panther_paint_dir

    main:

    /* XREFS:
    Aggregate matches across all members for each sequence --> single JSON with all matches for the batch
    Add signature and entry desc and names
    Add PAINT annotations (if panther is enabled)
    Add go terms (if enabled)
    Add pathways (if enabled)
    */
    XREFS(
        match_results,
        applications,
        db_releases,
        add_goterms,
        add_pathways,
        panther_paint_dir
    )

    REPRESENTATIVE_LOCATIONS(XREFS.out)
    // Collect all JSON files into a single channel so we don't have cocurrent writing to the output files
    ch_results = REPRESENTATIVE_LOCATIONS.out
        .map { meta, json -> json }
        .collect()

    emit:
    ch_results
}