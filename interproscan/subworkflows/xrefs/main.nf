include { ENTRIES } from "$projectDir/interproscan/modules/xrefs/entries/main"
include { GOTERMS } from "$projectDir/interproscan/modules/xrefs/goterms/main"
include { PAINT_ANNOTATIONS } from "$projectDir/interproscan/modules/xrefs/paint_annotations/main"
include { PATHWAYS } from "$projectDir/interproscan/modules/xrefs/pathways/main"

workflow XREFS {
    take:
    matches
    applications

    main:
    log.info "${params.xrefs.entries}"
    ENTRIES(matches, params.xrefs.entries)

    final_result = ENTRIES.out

    if ("${applications}".contains('panther')) {
        PAINT_ANNOTATIONS(final_result, params.members.panther.postprocess.paint_annotations)
        final_result = PAINT_ANNOTATIONS.out
    }

    if (params.goterms) {
        GOTERMS(final_result, params.xrefs.goterms)
        final_result = GOTERMS.out
    }

    if (params.pathways) {
        PATHWAYS(final_result, params.xrefs.pathways)
        final_result = PATHWAYS.out
    }

    emit:
    final_result
}
