import groovy.json.JsonOutput
import groovy.json.JsonSlurper

include {
    ENTRIES;
    GOTERMS;
    PATHWAYS;
    PAINT_ANNOTATIONS;
    AGGREGATE_XREFS;
} from "$projectDir/interproscan/modules/xrefs/main"


workflow XREFS {
    take:
    matches
    dataDir
    applications

    main:
    def entries_path = "${dataDir}/${params.xrefs.entries}"
    ENTRIES(entries_path, matches)

    matches_paint = Channel.empty()
    goterms_output = Channel.empty()
    pathways_output = Channel.empty()

    if ("${applications}".contains('panther')) {
        def paint_anno_dir = "${dataDir}/${params.members."panther".postprocess.paint_annotations}"
        matches_paint = PAINT_ANNOTATIONS(paint_anno_dir, matches2entries)
    }

    if (params.goterms || params.pathways) {
        matches2interpro = ENTRIES.out.matchkey2interpro

        if (params.goterms) {
            def ipr2go_path = "${dataDir}/${params.xrefs.goterms}.ipr.json"
            def go_info_path = "${dataDir}/${params.xrefs.goterms}.json"
            goterms_output = GOTERMS(ipr2go_path, go_info_path, matches2interpro)
        }
        if (params.pathways) {
            def ipr2pa_path = "${dataDir}/${params.xrefs.pathways}.ipr.json"
            def pa_info_path = "${dataDir}/${params.xrefs.pathways}.json"
            pathways_output = PATHWAYS(ipr2pa_path, pa_info_path, matches2interpro)
        }
    }

    AGGREGATE_XREFS(
        goterms_output.collect(),
        pathways_output.collect(),
        ENTRIES.out.matches2entries.collect(),
        matches_paint.collect()
    )

    emit:
    AGGREGATE_XREFS.out
}
