import groovy.json.JsonOutput
import groovy.json.JsonSlurper

process PAINT_ANNOTATIONS {
    label 'xref'

    // Retrieve PAINT annotations for Panther hits
    // calculated and pre-calc becuase they are not retrieved from the Match Lookup
    input:
    val paint_anno_dir
    val ch_matches2xrefs

    output:
    path matches_paintanno

    exec:
    def outputFilePath = task.workDir.resolve("matches_paintanno.json")
    def matches_annot = [:]
    matches_annot = ch_matches2xrefs.collectEntries { seq_id, match_info ->
        match_info.each { sig_acc, data ->
            if (data["member_db"].toUpperCase() == "PANTHER") {
                def anno_path = "${paint_anno_dir}/${sig_acc}.json"
                def paint_annotation_file = new File(anno_path)
                if (paint_annotation_file.exists()) {
                    def paint_annotations_content = new JsonSlurper().parse(paint_annotation_file)
                    def node_data = paint_annotations_content[data["node_id"]]
                    data["proteinClass"] = node_data[2]
                    data["graftPoint"] = node_data[3]
                }
            }
        }
        return [(seq_id): match_info]
    }
    def json = JsonOutput.toJson(matches_annot)
    new File(outputFilePath.toString()).write(json)
}
