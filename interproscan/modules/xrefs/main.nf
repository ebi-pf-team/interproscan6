import groovy.json.JsonOutput
import groovy.json.JsonSlurper

process GOTERMS {
    label 'xref'

    input:
    val ipr2go_path
    val go_info_path
    val matchesinterpro

    output:
    path "match_key2go.txt"

    script:
    def GO_PATTERN = [
        "P": "BIOLOGICAL_PROCESS",
        "C": "CELLULAR_COMPONENT",
        "F": "MOLECULAR_FUNCTION"
    ]

    def ipr2go = new JsonSlurper().parseText(new File(ipr2go_path).text)
    def go_info = new JsonSlurper().parseText(new File(go_info_path).text)
    def matches_info = matchesinterpro[1]
    def match_key2go = [:]
    matches_info.each { match_key, interpro_key ->
        try {
            def go_ids = ipr2go[interpro_key]
            for (String go_id : go_ids) {
                def go_dict = [
                    "name": go_info[go_id][1],
                    "databaseName": GO_PATTERN[go_info[go_id][0]],
                    "id": go_id
                ]
                if (!match_key2go.containsKey(match_key)) {
                    match_key2go[match_key] = []
                }
                match_key2go[match_key] << go_dict
            }
        } catch (Exception e) {
            match_key2go[match_key] = []
        }
    }
    """
    echo "${match_key2go}" > match_key2go.txt
    """
}

process PATHWAYS {
    label 'xref'

    input:
    val ipr2pa_path
    val pa_info_path
    val matchesinterpro

    output:
    path "match_key2pa.txt"

    script:
    def PA_PATTERN = [
    "t": "MetaCyc",
    "w": "UniPathway",
    "k": "KEGG",
    "r": "Reactome"
    ]

    def ipr2pa = new JsonSlurper().parseText(new File(ipr2pa_path).text)
    def pa_info = new JsonSlurper().parseText(new File(pa_info_path).text)
    def matches_info = matchesinterpro[1]
    def match_key2pa = [:]
    matches_info.each { match_key, interpro_key ->
        try {
            def pa_ids = ipr2pa[interpro_key]
            for (String pa_id : pa_ids) {
                def pa_dict = [
                    "name": pa_info[pa_id][1],
                    "databaseName": PA_PATTERN[pa_info[pa_id][0]],
                    "id": pa_id
                ]
                if (!match_key2pa.containsKey(match_key)) {
                    match_key2pa[match_key] = []
                }
                match_key2pa[match_key] << pa_dict
            }
        } catch (Exception e) {
            match_key2pa[match_key] = []
        }
    }
    """
    echo "${match_key2pa}" > match_key2pa.txt
    """
}

// process PAINT_ANNOTATIONS {
//     label 'xref'
//
//     input:
//     val paint_anno_dir
//     val ch_matches2xrefs
//
//     output:
//     val paint_annotations
//
//     script:
//     ch_matches2xrefs = ch_matches2xrefs.collectEntries { seq_id, match_info ->
//         match_info.each { sig_acc, data ->
//             if (data["member_db"].toUpperCase() == "PANTHER") {
//                 def anno_path = "${paint_anno_dir}/${sig_acc}.json"
//                 def paint_annotation_file = new File(anno_path)
//                 if (paint_annotation_file.exists()) {
//                     def paint_annotations_content = new JsonSlurper().parse(paint_annotation_file)
//                     def node_data = paint_annotations_content[data["node_id"]]
//                     data["proteinClass"] = node_data[2]
//                     data["graftPoint"] = node_data[3]
//             }
//         }
//         return [(seq_id): match_info]
//     }
// }
