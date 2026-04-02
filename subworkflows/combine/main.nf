include { COMBINE_MATCHES; COMBINE_MATCHES_BULK                   } from "../../modules/combine"
include { XREFS; XREFS_BULK                                       } from "../../modules/xrefs"
include { REPRESENTATIVE_LOCATIONS; REPRESENTATIVE_LOCATIONS_BULK } from "../../modules/representative_locations"

workflow COMBINE {
    take:
    json                // [ meta, [json] ]
    appl_config
    appl_dirs
    add_goterms
    add_pathways
    skip_repr_locations
    batch_size

    main:
    ch_interpro = appl_dirs
        .filter { name, dirpath -> name == "interpro" }
        .map    { name, dirpath -> dirpath }
        .ifEmpty(file("${projectDir}/assets/DUMMY")) 
        .first()

    ch_panther = appl_dirs
        .filter { name, dirpath -> name == "panther" }
        .map    { name, dirpath -> dirpath.resolve(appl_config.panther.paint) }
        .ifEmpty(file("${projectDir}/assets/DUMMY2")) 
        .first()

    ch_xrefs = channel.empty()
    if (batch_size < 50000) {
        COMBINE_MATCHES(json)

        ch_xrefs = XREFS(
            COMBINE_MATCHES.out,
            ch_interpro,
            ch_panther,
            add_goterms,
            add_pathways
        )
    } else {
        COMBINE_MATCHES_BULK(json)

        ch_xrefs = XREFS_BULK(
            COMBINE_MATCHES_BULK.out,
            ch_interpro,
            ch_panther,
            add_goterms,
            add_pathways
        )
    }

    ch_combined = channel.empty()
    if (skip_repr_locations) {
        ch_combined = ch_xrefs
    } else if (batch_size < 50000) {
        ch_combined = REPRESENTATIVE_LOCATIONS(ch_xrefs)
    } else {
        ch_combined = REPRESENTATIVE_LOCATIONS_BULK(ch_xrefs)
    }

    // Collect all JSON files into a single channel so we don't have cocurrent writing to the output files
    results = ch_combined
        .map { meta, json -> json }
        .collect()

    emit:
    results
}