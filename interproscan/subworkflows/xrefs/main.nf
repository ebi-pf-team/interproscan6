import groovy.json.JsonOutput
import groovy.json.JsonSlurper

include {
    PANTHER_XREFS;
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
            }
        }
        return matches_info
    }

    if ("${applications}".contains('panther')) {
        def paint_annotations = new JsonSlurper().parseText(new File("${dataDir}/${params.members."panther".postprocess.paint_annotations}").text)
        matches2annotations = PANTHER_XREFS(paint_annotations, matches2entries)
    }

    if (params.goterms) {
        def ipr2go = new JsonSlurper().parseText(new File("${dataDir}/${params.xrefs.goterms}.ipr.json").text)
        def go_info = new JsonSlurper().parseText(new File("${dataDir}/${params.xrefs.goterms}.json").text)
        GOTERMS(ipr2go, go_info, matches2entries)
    }

    if (params.pathways) {
        def ipr2pa = new JsonSlurper().parseText(new File("${dataDir}/${params.xrefs.pathways}.ipr.json").text)
        def pa_info = new JsonSlurper().parseText(new File("${dataDir}/${params.xrefs.pathways}.json").text)
        PATHWAYS(ipr2pa, pa_info, matches2entries)
    }

    ch_aggregated_results = matches2entries
        .collect()
        .map { jsonObjects ->
            def merged = [:]
            jsonObjects.each { hits ->
                hits.each { key, value ->
                    merged[key] = merged.containsKey(key) ? merged[key] + value : value
                }
            }
            def aggregated_result = JsonOutput.toJson(merged)
            def outputFile = new File("${workDir}/aggregated_result.json")
            outputFile.write(aggregated_result)
            return outputFile.path
        }

    emit:
    ch_aggregated_results
}
