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
            }
        }
        return [matches_info, matches2interpro]
    }

//     if ("${applications}".contains('panther')) {
//         def paint_annotations = new JsonSlurper().parseText(new File("${dataDir}/${params.members."panther".postprocess.paint_annotations}").text)
//         PANTHER_XREFS(paint_annotations, matches2entries).set { panther_output }
//     }

    if (params.goterms || params.pathways) {
        def goterms_output = Channel.value([[:]])
        def pathways_output = Channel.value([[:]])

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

        ch_combined = matches2entries[0].combine(goterms_output, pathways_output)
            .map { matches, goXRefs, pathwayXRefs ->
                def matches2xrefs = [:]
                matches.each { hits ->
                    hits.each { match_key, value ->
                        matches2xrefs[match_key] = matches2xrefs.containsKey(match_key) ? matches2xrefs[match_key] + value : value
                        matches2xrefs[match_key]["entry"]["goXRefs"] = goXRefs[match_key] ?: []
                        matches2xrefs[match_key]["entry"]["pathwayXRefs"] = pathwayXRefs[match_key] ?: []
                    }
                }
            return matches2xrefs
        }
    } else {
        ch_combined = matches2entries[0]
    }

    ch_aggregated_results = ch_combined.collect().map { all_matches2xrefs ->
        def aggregated_result = JsonOutput.toJson(all_matches2xrefs)
        def outputFile = new File("${workDir}/aggregated_result.json")
        outputFile.write(aggregated_result)
        return outputFile.path
    }

    emit:
    ch_aggregated_results
}
