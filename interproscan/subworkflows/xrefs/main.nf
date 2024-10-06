import groovy.json.JsonOutput
import groovy.json.JsonSlurper

include {
    GOTERMS;
    PATHWAYS;
    PAINT_ANNOTATIONS;
    AGGREGATE_XREFS;
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

    matches_paint = Channel.empty()
    goterms_output = Channel.empty()
    pathways_output = Channel.empty()

    if ("${applications}".contains('panther')) {
        def paint_anno_dir = "${dataDir}/${params.members."panther".postprocess.paint_annotations}"
        matches_paint = PAINT_ANNOTATIONS(paint_anno_dir, matches2entries)
    }

    if (params.goterms || params.pathways) {
        if (params.goterms) {
            def ipr2go_path = "${dataDir}/${params.xrefs.goterms}.ipr.json"
            def go_info_path = "${dataDir}/${params.xrefs.goterms}.json"
            goterms_output = GOTERMS(ipr2go_path, go_info_path, matches2entries)
        }
        if (params.pathways) {
            def ipr2pa_path = "${dataDir}/${params.xrefs.pathways}.ipr.json"
            def pa_info_path = "${dataDir}/${params.xrefs.pathways}.json"
            pathways_output = PATHWAYS(ipr2pa_path, pa_info_path, matches2entries)
        }
    }

    AGGREGATE_XREFS(
        goterms_output.collect(),
        pathways_output.collect(),
        matches2entries.collect(),
        matches_paint.collect()
    )

    emit:
    AGGREGATE_XREFS.out
}
