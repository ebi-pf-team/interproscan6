import groovy.json.JsonOutput
import groovy.json.JsonSlurper

include {
    GOTERMS;
    PATHWAYS;
//     PAINT_ANNOTATIONS;
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
                println "match_key: ${match_key}"
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

//     if ("${applications}".contains('panther')) {
//         def paint_anno_dir = "${dataDir}/${params.members."panther".postprocess.paint_annotations}"
//         matches_paint = PAINT_ANNOTATIONS(paint_anno_dir, ch_matches2xrefs)
//     }

    if (params.goterms || params.pathways) {
        def goterms_output = Channel.empty()
        def pathways_output = Channel.empty()
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

        def go_pa_output = goterms_output.mix(pathways_output).collect()
        ch_matches2xrefs = matches2entries.map { entries_output ->
            def matches_info = entries_output[0]
            return matches_info.collectEntries { seq_id, match_data ->
                def match_key = match_data.keySet().first()
                def match_info = match_data[match_key]

                match_info_xrefs = go_pa_output.map { file ->
                    def xRefsMap = new File(file.toString()).text
                    if file.toString().contains("go") {
                         def goXRefs_info = goRefsMap[match_key] ?: []
                         match_info["entry"]["goXRefs"] = goXRefs_info
                    } else {
                         def paXRefs_info = paRefsMap[match_key] ?: []
                         match_info["entry"]["pathwayXRefs"] = paXRefs_info
                    }
                    return match_info
                }
                return [(seq_id): [(match_key): match_info_xrefs]])
            }
        }
    } else {
        ch_matches2xrefs = matches2entries.collect().map { entries_output ->
            def matches_info = entries_output[0]
            return matches_info
        }
    }

//     def outputFile = new File("${workDir}/aggregated_result.json")
//     outputFile.write(ch_matches2xrefs.collect())

    emit:
    ch_matches2xrefs
}
