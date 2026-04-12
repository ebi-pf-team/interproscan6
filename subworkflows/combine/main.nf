include { COMBINE_MATCHES          } from "../../modules/combine"
include { ADD_XREFS                } from "../../modules/add_xrefs"
include { REPRESENTATIVE_LOCATIONS } from "../../modules/representative_locations"

workflow COMBINE {
    take:
    json                // [ meta, [json] ]
    appl_config
    appl_dirs
    add_goterms
    add_pathways
    skip_repr_locations

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

    COMBINE_MATCHES(json)

    ch_xrefs = ADD_XREFS(
        COMBINE_MATCHES.out,
        ch_interpro,
        ch_panther,
        add_goterms,
        add_pathways
    )

    ch_combined = channel.empty()
    if (skip_repr_locations) {
        ch_combined = ch_xrefs
    } else {
        ch_combined = REPRESENTATIVE_LOCATIONS(ch_xrefs)
    }

    // Collect all JSON files into a single channel so we don't have concurrent writing to the output files
    results = ch_combined
        .map { meta, json -> json }
        .collect(flat: true)

    emit:
    results
}