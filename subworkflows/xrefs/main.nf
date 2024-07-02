include { ENTRIES } from "$projectDir/modules/xrefs/entries/main"
include { GOTERMS } from "$projectDir/modules/xrefs/goterms/main"
// include { PAINT_ANNOTATIONS } from "$projectDir/modules/xref/paint_annotations/main"
include { PATHWAYS } from "$projectDir/modules/xrefs/pathways/main"

workflow XREFS {
    take:
    matches
    applications

    main:

    ENTRIES(matches, params.xrefs.entries)

    // PAINT_ANNOTATIONS(matches, params.members.panther.postprocess.paint_annotations)

    if (params.goterms) {
        GOTERMS(ENTRIES.out, params.xrefs.goterms)
        if (params.pathways) {
            final_result = PATHWAYS(GOTERMS.out, params.xrefs.pathways)
        }
        else {
            final_result = GOTERMS.out
        }
    }
    else {
        if (params.pathways) {
            final_result = PATHWAYS(ENTRIES.out, params.xrefs.pathways)
        }
        else {
            final_result = ENTRIES.out
        }
    }

    emit:
    final_result
}
