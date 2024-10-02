include {
    GOTERMS;
    PATHWAYS;
} from "$projectDir/interproscan/modules/xrefs/main"

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

workflow XREFS {
    take:
    matches
    dataDir

    main:
    def entries = new JsonSlurper().parseText(new File("${dataDir}/${params.xrefs.entries}").text)
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
                    def paint_annotations = new JsonSlurper().parseText(new File("${dataDir}/${params.members."panther".postprocess.paint_annotations}/{sig_acc}.json").text)
                    def acc_subfamily = data["accession"]
                    try {
                        match_info[match_key]["entry"]["subfamily_name"] = entries_info[acc_subfamily]["name"]
                        match_info[match_key]["entry"]["subfamily_description"] = entries_info[acc_subfamily]["description"]
                        match_info[match_key]["entry"]["subfamily_type"] = entries_info[acc_subfamily]["type"]
                    } catch (Exception e) {
                        log.info "No subfamily info for ${acc_subfamily}"
                    }
                    node_data = paint_annotations[data["node_id"]]
                    match_info[seq_id][sig_acc]["proteinClass"] = node_data[2]
                    match_info[seq_id][sig_acc]["graftPoint"] = node_data[3]
                }
            }
        }
        return matches_info
    }

    if (params.goterms) {
        def ipr2go = new JsonSlurper().parseText(new File("${dataDir}/${params.xrefs.goterms}.ipr.json").text)
        def go_info = new JsonSlurper().parseText(new File("${dataDir}/${params.xrefs.goterms}.json").text)
        matches2entries = GOTERMS(ipr2go, go_info, matches2entries)
    }
    if (params.pathways) {
        def ipr2pa = new JsonSlurper().parseText(new File("${dataDir}/${params.xrefs.pathways}.ipr.json").text)
        def pa_info = new JsonSlurper().parseText(new File("${dataDir}/${params.xrefs.pathways}.json").text)
        matches2entries = PATHWAYS(ipr2pa, pa_info, matches2entries)
    }

    emit:
    matches2entries
}
