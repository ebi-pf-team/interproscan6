import groovy.json.JsonSlurper

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

def GO_PATTERN = [
    "P": "BIOLOGICAL_PROCESS",
    "C": "CELLULAR_COMPONENT",
    "F": "MOLECULAR_FUNCTION"
]

def PA_PATTERN = [
    "t": "MetaCyc",
    "w": "UniPathway",
    "k": "KEGG",
    "r": "Reactome"
]


workflow XREFS {
    take:
    matches
    dataDir

    main:
    def matches2entries = [:]
    def matches_info = new JsonSlurper().parseText(new File(matches_path).text)
    def entries = new JsonSlurper().parseText(new File("${dataDir}/${params.xrefs.entries}").text)

    matches_info.each { seq_id, match_info ->
        match_info.each { match_key, data ->
            def databases_versions = entries["databases"]
            def entries_info = entries['entries']
            def acc_id = match_key.split("\\.")[0]
            def member_db = data["member_db"].toUpperCase()
            match_info[match_key]['member_db'] = member_db
            match_info[match_key]['library'] = MAP_DATABASES[member_db]
            if (!match_info[match_key].containsKey('version')) {  // mobidb, signalp and tmhmm get version from members.config
                match_info[match_key]['version'] = databases_versions[map_databases[member_db]]
            }

            def entry = entries_info[acc_id] ?: entries_info[match_key]
            if (entry) {
                match_info[match_key]['entry'] = entry
                def interpro_key = entry['integrated']
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
            }
            else {  // members with no match in entries. e.g. PIRSR
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
                    log.info "No subfamily info for ${acc_subfamily}"
                }
                def paint_annotations = new JsonSlurper().parseText(new File("${dataDir}/${params.members."panther".postprocess.paint_annotations}/{sig_acc}.json").text)
                node_data = paint_annotations[data["node_id"]]
                matches[seq_id][sig_acc]["proteinClass"] = node_data[2]
                matches[seq_id][sig_acc]["graftPoint"] = node_data[3]
                }
            }

            if (data["entry"]) {
                def ipr_id = data["entry"]["accession"]
                if (params.goterms) {
                    def ipr2go = new JsonSlurper().parseText(new File("${dataDir}/${params.xrefs.goterms.ipr.json}").text)
                    def go_info = new JsonSlurper().parseText(new File("${dataDir}/${params.xrefs.goterms.json}").text)
                    try {
                        go_ids = ipr2go[ipr_id]
                        for (String go_id : go_ids) {
                            go_dict = [
                                "name": go_info[go_id][0],
                                "databaseName": "GO",
                                "category": db_pattern[go_info[go_id][1]],
                                "id": go_id
                            ]
                            match_info[match_key]["entry"]["goXRefs"] + go_dict
                        }
                    } catch (Exception e) {
                        log.info "No GO terms for ${ipr_id}"
                    }

                if (params.pathways) {
                    def ipr2pa = new JsonSlurper().parseText(new File("${dataDir}/${params.xrefs.pathways.ipr.json}").text)
                    def pa_info = new JsonSlurper().parseText(new File("${dataDir}/${params.xrefs.pathways.json}").text)
                    try {
                        pa_ids = ipr2pa[ipr_id]
                        for (String pa_id : pa_ids) {
                            pa_dict = [
                                "name": pa_info[pa_id][1],
                                "databaseName": db_pattern[pa_info[pa_id][0]],
                                "id": pa_id
                            ]
                            match_info[match_key]["entry"]["pathwayXRefs"] + pa_dict
                        }
                    } catch (Exception e) {
                        log.info "No Pathways for ${ipr_id}"
                    }
                }
            }
        }
        matches2entries[seq_id] = match_info
    }

    emit:
    matches2entries
}
