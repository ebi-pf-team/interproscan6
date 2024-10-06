import groovy.json.JsonOutput
import groovy.json.JsonSlurper

process GOTERMS {
    label 'xref'

    input:
    val ipr2go_path
    val go_info_path
    val matchesinterpro

    output:
    path "match_key2go.json"

    exec:
    def GO_PATTERN = [
        "P": "BIOLOGICAL_PROCESS",
        "C": "CELLULAR_COMPONENT",
        "F": "MOLECULAR_FUNCTION"
    ]

    def outputFilePath = task.workDir.resolve("match_key2go.json")
    def ipr2go = new JsonSlurper().parseText(new File(ipr2go_path).text)
    def go_info = new JsonSlurper().parseText(new File(go_info_path).text)
    def matches_info = matchesinterpro[1]
    def match_key2go = [:]
    match_key2go = matches_info.collectEntries { match_key, interpro_key ->
        try {
            def go_ids = ipr2go[interpro_key]
            def go_terms = go_ids.collect { go_id ->
                [
                    "name": go_info[go_id][1],
                    "databaseName": GO_PATTERN[go_info[go_id][0]],
                    "id": go_id
                ]
            }
            return [(match_key): go_terms]
        } catch (Exception e) {
            return [(match_key): []]
        }
    }
    def json = JsonOutput.toJson(match_key2go)
    new File(outputFilePath.toString()).write(json)
}

process PATHWAYS {
    label 'xref'

    input:
    val ipr2pa_path
    val pa_info_path
    val matchesinterpro

    output:
    path "match_key2pa.json"

    exec:
    def PA_PATTERN = [
    "t": "MetaCyc",
    "w": "UniPathway",
    "k": "KEGG",
    "r": "Reactome"
    ]

    def outputFilePath = task.workDir.resolve("match_key2pa.json")
    def ipr2pa = new JsonSlurper().parseText(new File(ipr2pa_path).text)
    def pa_info = new JsonSlurper().parseText(new File(pa_info_path).text)
    def matches_info = matchesinterpro[1]

    def match_key2pa = [:]
    match_key2pa = matches_info.collectEntries { match_key, interpro_key ->
        try {
            def pa_ids = ipr2pa[interpro_key]
            def pa_terms = pa_ids.collect { pa_id ->
                [
                    "name": pa_info[pa_id][1],
                    "databaseName": PA_PATTERN[pa_info[pa_id][0]],
                    "id": pa_id
                ]
            }
            return [(match_key): pa_terms]
        } catch (Exception e) {
            return [(match_key): []]
        }
    }
    def json = JsonOutput.toJson(match_key2pa)
    new File(outputFilePath.toString()).write(json)
}

process PAINT_ANNOTATIONS {
    label 'xref'

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

process AGGREGATE_XREFS {
    label 'xref'

    input:
    val goterms_output
    val pathways_output
    val matches2entries
    val paint_annotations_output

    output:
    path "aggregated_matches_xrefs.json"

    exec:
    def go_pa_output = goterms_output.mix(pathways_output)
    println "GO PA OUTPUT: ${go_pa_output}"
    def aggregated_matches_xrefs = [:]
    aggregated_matches_xrefs = matches2entries.map { entries_output ->
        def matches_info = entries_output[0]
        return matches_info.collectEntries { seq_id, match_data ->
            def match_key = match_data.keySet().first()
            def match_info = match_data[match_key]
            def match_info_xrefs = go_pa_output.collectEntries { file ->
                println "File: ${file}"
                print "File: ${file.getName()}"
                if (file.getName().contains('match_key2go')) {
                    println(file)
//                     match_info["entry"]["goXRefs"] = "match_key2go"
//                     match_info["entry"]["goXRefs"] = goterms ?: []
                } else {
                    println(file)
//                     match_info["entry"]["pathwayXRefs"] = "match_key2pa"
//                     match_info["entry"]["pathwayXRefs"] = pathways ?: []
                }
                return match_info
            }
            return [(seq_id): [(match_key): match_info_xrefs]]
        }
    }
    new File("aggregated_matches_xrefs.json").write(JsonOutput.toJson(aggregated_matches_xrefs))
}
