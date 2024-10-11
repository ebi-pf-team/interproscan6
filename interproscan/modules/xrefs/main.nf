import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import com.fasterxml.jackson.databind.ObjectMapper

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

process XREFS {
    label 'xrefs'

    input:
    tuple val(meta), val(membersMatches)
    val apps
    val data_dir

    output:
    tuple val(meta), path("matches2xrefs.json")

    exec:
    // Getting xrefs data files
    String entriesPath = "${data_dir}/${params.xrefs.entries}"
    File entriesJson = new File(entriesPath.toString())
    def objectMapper = new ObjectMapper()
    def entries = objectMapper.readValue(entriesJson, Map)

    JsonSlurper jsonSlurper = new JsonSlurper()

    if ("${apps}".contains('panther')) {
        String paintAnnoDir = "${data_dir}/${params.xrefs.paint_annotations}"
    }

    if (params.goterms) {
        String ipr2goPath = "${data_dir}/${params.xrefs.goterms}.ipr.json"
        String goInfoPath = "${data_dir}/${params.xrefs.goterms}.json"
        def ipr2go = jsonSlurper.parse(new File(ipr2goPath.toString()).text)
        def goInfo = jsonSlurper.parse(new File(goInfoPath.toString()).text)
    }

    if (params.pathways) {
        String ipr2paPath = "${data_dir}/${params.xrefs.pathways}.ipr.json"
        String paInfoPath = "${data_dir}/${params.xrefs.pathways}.json"
        def ipr2pa = jsonSlurper.parseText(new File(ipr2paPath.toString()).text)
        def paInfo = jsonSlurper.parseText(new File(paInfoPath.toString()).text)
    }

    // Enriching matches with xrefs
    matchesEntries = membersMatches.each { matchesPath  ->
        memberDB = matchesPath.toString().split("/").last().split("\\.")[0]
        def sigLibRelease = [
            "library": memberDB,
            "version": entries['databases'][memberDB]
        ]
        def matches = jsonSlurper.parse(matchesPath).collectEntries { seqId, jsonMatches ->
            [(seqId): jsonMatches.collectEntries { matchId, jsonMatch ->
                Match matchObject = Match.fromMap(jsonMatch)
                def entriesInfo = entries['entries']
                String accId = matchObject.modelAccession
                def entry = entriesInfo[accId] ?: entriesInfo[matchKey]
                entryData = null
                if (entry) {
                    def interproKey = entry['integrated']
                    def entryInfo = entriesInfo.get(interproKey)
                    if (entryInfo) {
                        entryData = [
                            "accession": interproKey,
                            "name": entryInfo["name"],
                            "description": entryInfo["description"],
                            "type": entryInfo["type"]
                        ]
                    }
                    def representativeFlag = null
                    if (entry['representative']) {
                        representative = [
                            "type": entry['representative']["type"],
                            "rank": entry['representative']["rank"]
                        ]
                        RepresentativeInfo representativeInfo = RepresentativeInfo.fromMap(representative)
                    }
                    matchObject.representativeInfo = representativeInfo

                    if (matchObject.signature.signatureLibraryRelease.library == "panther") {
                        String sigAcc = matchObject.signature.accession
                        String paintAnnPath = "${paintAnnDir}/${sigAcc}.json"
                        File paintAnnotationFile = new File(paintAnnPath.toString())
                        if (paintAnnotationFile.exists()) {
                            def paintAnnotationsContent = jsonSlurper.parse(paintAnnotationFile)
                            String nodeId = matchObject.treegrafter.ancestralNodeID
                            def nodeData = paintAnnotationsContent[nodeId]
                            matchObject.treegrafter.proteinClass = nodeData[2]
                            matchObject.treegrafter.graftPoint = nodeData[3]
                        }
                    }

                    if (params.goterms || params.pathways) {
                        String interproKey = matchObject.signature.entry.accession
                        try {
                            if (params.goterms) {
                                def goIds = ipr2go[interproKey]
                                def goTerms = goIds.collect { goId ->
                                    goXref = [
                                        "name": goInfo[goId][0],
                                        "databaseName": "GO",
                                        "category": GO_PATTERN[goInfo[goId][1]],
                                        "id": goId
                                    ]
                                    GoXrefs goXrefsObj = GoXrefs.fromMap(goXref)
                                    matchObject.signature.entry.goXrefs.add(goXrefsObj)
                                }
                            }
                            if (params.pathways) {
                                def paIds = ipr2pa[interproKey]
                                def paTerms = paIds.collect { paId ->
                                    paXref = [
                                        "name": paInfo[paId][1],
                                        "databaseName": PA_PATTERN[paInfo[paId][0]],
                                        "id": paId
                                    ]
                                    PathwayXrefs pathwayXrefsObj = PathwayXrefs.fromMap(paXref)
                                    matchObject.signature.entry.pathwayXrefs.add(pathwayXrefsObj)
                                }
                            }
                        } catch (Exception e) {
                            // pass
                        }
                    }
                }

                def sigData = [
                    "accession": entry["accession"],
                    "name": entry["name"],
                    "description": entry["description"],
                    "signatureLibraryRelease": sigLibRelease,
                    "entry": entryData
                ]
                Signature signatureObject = Signature.fromMap(sigData)
                matchObject.signature = signatureObject

                if (memberDB == "panther") {
                    String accSubfamily = matchObject.treegrafter.subfamilyAccession
                    matchObject.treegrafter.subfamilyAccession = accSubfamily
                    if (entriesInfo[accSubfamily]) {
                        matchObject.treegrafter.subfamilyName = entrieInfo[accSubfamily]["name"]
//                         Just waiting finish output step to have sure we don't use this infos
//                         matchObject.treegrafter.subfamilyDesc = entriesInfo[accSubfamily]["description"]
//                         matchObject.treegrafter.subfamilyType = entriesInfo[accSubfamily]["type"]
                    }
                }

                [(matchId): matchObject]
            }]
        }
        def outputFilePath = task.workDir.resolve("matches2xrefs.json")
        def json = JsonOutput.toJson(matches)
        new File(outputFilePath.toString()).write(json)
    }
}
