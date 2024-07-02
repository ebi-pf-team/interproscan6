include { ENTRIES } from "$projectDir/modules/xrefs/entries/main"
include { GOTERMS } from "$projectDir/modules/xrefs/goterms/main"
include { PAINT_ANNOTATIONS } from "$projectDir/modules/xrefs/paint_annotations/main"
include { PATHWAYS } from "$projectDir/modules/xrefs/pathways/main"

workflow XREFS {
    take:
    matches
    applications

    main:

    ENTRIES(matches, params.xrefs.entries)

    if ("${applications}".contains('panther')) {
        PAINT_ANNOTATIONS(ENTRIES.out, params.members.panther.postprocess.paint_annotations)
    }

    if (params.goterms) {
        if ("${applications}".contains('panther')) {
            GOTERMS(PAINT_ANNOTATIONS.out, params.xrefs.goterms)
        }
        else {
            GOTERMS(ENTRIES.out, params.xrefs.goterms)
        }
        
        if (params.pathways) {
            final_result = PATHWAYS(GOTERMS.out, params.xrefs.pathways)
        }
        else {
            final_result = GOTERMS.out
        }
    }
    else {
        if (params.pathways) {
            if ("${applications}".contains('panther')) {
                final_result = PATHWAYS(PAINT_ANNOTATIONS.out, params.xrefs.pathways)
            }
            else {
                final_result = PATHWAYS(ENTRIES.out, params.xrefs.pathways)
            }
        }
        else {
            if ("${applications}".contains('panther')) {
                final_result = PAINT_ANNOTATIONS.out
            }
            else {
                final_result = ENTRIES.out
            }
        }
    }

    emit:
    final_result
}
