import groovy.json.JsonSlurper

include { ENTRIES } from "../../modules/xrefs/entries"
include { GOTERMS } from "../../modules/xrefs/goterms"
include { PAINT_ANNOTATIONS } from "../../modules/xrefs/paint_annotations"
include { PATHWAYS } from "../../modules/xrefs/pathways"
include { AGGREGATE_RESULTS } from "../../modules/xrefs/aggregate_results"

workflow XREFS {
    take:
    matches
    apps
    data_dir

    main:
    def entriesPath = "${data_dir}/${params.xrefs.entries}"
    File entriesJson = new File(entriesPath.toString())
    ENTRIES(matches, entriesJson)

    matches_paint = Channel.empty()
    matchesGoObj = Channel.empty()
    matchesPaObj = Channel.empty()

    if ("${apps}".contains('panther')) {
        def paint_anno_dir = "${data_dir}/${params.members."panther".postprocess.paint_annotations}"
        matches_paint = PAINT_ANNOTATIONS(paint_anno_dir, ENTRIES.out)
    }

    if (params.goterms) {
        def ipr2goPath = "${data_dir}/${params.xrefs.goterms}.ipr.json"
        def goInfoPath = "${data_dir}/${params.xrefs.goterms}.json"
        File ipr2goJson = new File(ipr2goPath.toString())
        File goInfoJson = new File(goInfoPath.toString())
        matchesGoObj = GOTERMS(ipr2goJson, goInfoJson, ENTRIES.out)
    }

    if (params.pathways) {
        def ipr2paPath = "${data_dir}/${params.xrefs.pathways}.ipr.json"
        def paInfoPath = "${data_dir}/${params.xrefs.pathways}.json"
        File ipr2paJson = new File(ipr2paPath.toString())
        File paInfoJson = new File(paInfoPath.toString())
        matchesPaObj = PATHWAYS(ipr2paJson, paInfoJson, ENTRIES.out)
    }

    AGGREGATE_RESULTS(ENTRIES.out, matches_paint, matchesGoObj, matchesPaObj)

    emit:
    AGGREGATE_RESULTS.out
}
