import groovy.json.JsonOutput
import groovy.json.JsonSlurper

include {
    GOTERMS;
    PATHWAYS;
} from "$projectDir/interproscan/modules/xrefs/main"

def MAP_DATABASES = [
    "SFLD": "SFLD",
    "PRINTS": "PRINTS",
    "PFAM": "Pfam",
    "INTERPRO": "InterPro",
    "CDD": "CDD",
    "PROSITE_PROFILES": "PROSITE profiles",
    "NCBIFAM": "NCBIfam",
    "PROSITE_PATTERNS": "PROSITE patterns",
    "HAMAP": "HAMAP",
    "SMART": "SMART",
    "PIRSF": "PIRSF",
    "PANTHER": "PANTHER",
    "GENE3D": "CATH-Gene3D",
    "SUPERFAMILY": "SUPERFAMILY",
    "ANTIFAM": "AntiFam",
    "FUNFAM": "FunFam",
    "MOBIDB_LITE": "MobiDB Lite",
    "PHOBIUS": "Phobius",
    "SIGNALP": "SignalP",
    "SIGNALP_EUK": "SignalP_Euk",
    "TMHMM": "TMHMM",
    "COILS": "COILS",
    "PIRSR": "PIRSR"
]

workflow XREFS {
    take:
    matches
    dataDir
    applications

    main:
    def entries = new JsonSlurper().parseText(new File("${dataDir}/${params.xrefs.entries}").text)
    def matches2interpro = [:]
    matches2entries = matches.map { matchFile ->
        def matches_info = new JsonSlurper().parse(matchFile)
        matches_info.each { seq_id, match_info ->
            match_info.each { match_key, data ->
                def databases_versions = entries["databases"]
                def entries_info = entries['entries']
                def acc_id = match_key.split("\\.")[0]
                def member_db = data["member_db"].toUpperCase()
                match_info[match_key]['member_db'] = member_db
                match_info[match_key]['library'] = MAP_DATABASES[member_db]
                if (!match_info[match_key].containsKey('version')) {  // just mobidb, signalp and tmhmm get version from members.config
                    match_info[match_key]['version'] = databases_versions[MAP_DATABASES[member_db]]
                }

                def entry = entries_info[acc_id] ?: entries_info[match_key]
                if (entry) {
                    match_info[match_key]['entry'] = entry
                    def interpro_key = entry['integrated']
                    matches2interpro[match_key] = interpro_key
                    match_info[match_key]["entry"] = [
                        "accession": interpro_key,
                        "name": entry["name"],
                        "description": entry["description"],
                        "type": entry["type"],
                        "representative": entry["representative"],
                        "member_db": member_db,
                        "goXRefs": [],
                        "pathwayXRefs": []
                    ]
                    def ipr_info = entries_info.get(interpro_key)
                    if (ipr_info) {
                        match_info[match_key]["entry"]["ipr_name"] = ipr_info["name"]
                        match_info[match_key]["entry"]["ipr_description"] = ipr_info["description"]
                        match_info[match_key]["entry"]["ipr_type"] = ipr_info["type"]
                    }
                } else {  // members with no match in entries. e.g. PIRSR
                    match_info[match_key]["entry"] = [
                        "accession": None,
                        "name": None,
                        "description": None,
                        "type": None,
                        "database": member_db,
                        "goXRefs": [],
                        "pathwayXRefs": []
                    ]
                }

                if (member_db == "PANTHER") {
                    def acc_subfamily = data["accession"]
                    try {
                        match_info[match_key]["entry"]["subfamily_name"] = entries_info[acc_subfamily]["name"]
                        match_info[match_key]["entry"]["subfamily_description"] = entries_info[acc_subfamily]["description"]
                        match_info[match_key]["entry"]["subfamily_type"] = entries_info[acc_subfamily]["type"]
                    } catch (Exception e) {
                        log.info "No subfamily for ${acc_subfamily}"
                    }
                }
            }
        }
        return [matches_info, matches2interpro]
    }

    if (params.goterms || params.pathways) {
        def goterms_output = Channel.empty()
        def pathways_output = Channel.empty()
        if (params.goterms) {
            def ipr2go = new JsonSlurper().parseText(new File("${dataDir}/${params.xrefs.goterms}.ipr.json").text)
            def go_info = new JsonSlurper().parseText(new File("${dataDir}/${params.xrefs.goterms}.json").text)
            goterms_output = GOTERMS(ipr2go, go_info, matches2entries)
        }
        if (params.pathways) {
            def ipr2pa = new JsonSlurper().parseText(new File("${dataDir}/${params.xrefs.pathways}.ipr.json").text)
            def pa_info = new JsonSlurper().parseText(new File("${dataDir}/${params.xrefs.pathways}.json").text)
            pathways_output = PATHWAYS(ipr2pa, pa_info, matches2entries)
        }

        ch_matches2xrefs = matches2entries
        .map { entries_output ->
            def matches_info = entries_output[0]
            return matches_info.collectEntries { seq_id, match_data ->
                def match_key = match_data.keySet().first()
                def match_info = match_data[match_key]

                def goXRefs = goterms_output.map { it[match_key] ?: [] }
                def pathwayXRefs = pathways_output.map { it[match_key] ?: [] }

                match_info["entry"]["goXRefs"] = goXRefs
                match_info["entry"]["pathwayXRefs"] = pathwayXRefs

                return [(seq_id): [(match_key): match_info]]
            }
        }
    } else {
        ch_matches2xrefs = matches2entries.map { entries_output ->
            def matches_info = entries_output[0]
            return matches_info
        }
    }

    if ("${applications}".contains('panther')) {
        def paint_annotations = new JsonSlurper().parse(new File("${dataDir}/${params.members."panther".postprocess.paint_annotations}"))
        ch_matches2xrefs = ch_matches2xrefs.collectEntries { seq_id, match_info ->
            match_info.each { sig_acc, data ->
                if (data["member_db"].toUpperCase() == "PANTHER") {
                    def anno_path = "${paint_anno_dir}/${sig_acc}.json"
                    def paint_annotation_file = new File(anno_path)
                    if (paint_annotation_file.exists()) {
                        def paint_annotations_content = new JsonSlurper().parse(paint_annotation_file)
                        def node_data = paint_annotations_content[data["node_id"]]
                        data["proteinClass"] = node_data[2]
                        data["graftPoint"] = node_data[3]
                    } else {
                        log.info "No Panther annotation to ${sig_acc}"
                    }
                }
            }
            return [(seq_id): match_info]
        }
    }

    ch_aggregated_results = ch_matches2xrefs.collect().map { all_matches2xrefs ->
        def aggregated_result = JsonOutput.toJson(all_matches2xrefs)
        def outputFile = new File("${workDir}/aggregated_result.json")
        outputFile.write(aggregated_result)
        return outputFile.path
    }

    emit:
    ch_aggregated_results
}
