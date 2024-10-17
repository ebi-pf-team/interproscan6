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
    String entriesPath = "${data_dir}/${params.xRefsConfig.entries}"
    File entriesJson = new File(entriesPath.toString())
    def entries = new ObjectMapper().readValue(entriesJson, Map)
    JsonSlurper jsonSlurper = new JsonSlurper()

    matchesEntries = membersMatches.each { matchesPath  ->
        memberDB = matchesPath.toString().split("/").last().split("\\.")[0]
        def sigLibRelease = [
            "library": memberDB,
            "version": entries['databases'][memberDB]
        ]
        def matches = jsonSlurper.parse(matchesPath).collectEntries { seqId, jsonMatches ->
            [(seqId): jsonMatches.collectEntries { matchId, jsonMatch ->
                Match matchObject = Match.fromMap(jsonMatch)

                if (memberDB == "panther") {
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

                def entriesInfo = entries['entries']
                String accId = matchObject.modelAccession.split("\\.")[0]
                def sigData = [
                    "accession": accId,
                    "name": "",
                    "description": "",
                    "signatureLibraryRelease": sigLibRelease,
                    "entry": null
                ]
                Signature signatureObject = Signature.fromMap(sigData)
                matchObject.signature = signatureObject

                def entry = entriesInfo[accId] ?: entriesInfo[matchId]
                if (entry) {
                    matchObject.signature.name = entry["name"]
                    matchObject.signature.description = entry["description"]

                    def representativeFlag = null
                    if (entry['representative']) {
                        representative = [
                            "type": entry['representative']["type"],
                            "rank": entry['representative']["index"]
                        ]
                        RepresentativeInfo representativeInfo = RepresentativeInfo.fromMap(representative)
                        matchObject.representativeInfo = representativeInfo
                    }

                    def interproKey = entry['integrated']
                    def entryInfo = entriesInfo.get(interproKey)
                    if (entryInfo) {
                        entryData = [
                            "accession": interproKey,
                            "name": entryInfo["name"],
                            "description": entryInfo["description"],
                            "type": entryInfo["type"],
                            "goXRefs": [],
                            "pathwayXRefs": []
                        ]
                        Entry entryDataObj = Entry.fromMap(entryData)
                        matchObject.signature.entry = entryDataObj
                    }

                    if (params.goterms || params.pathways) {
                        interproKey = entry["integrated"]
                        if (params.goterms) {
                            try {
                                String ipr2goPath = "${data_dir}/${params.xRefsConfig.goterms}.ipr.json"
                                String goInfoPath = "${data_dir}/${params.xRefsConfig.goterms}.json"
                                ipr2go = jsonSlurper.parse(new File(ipr2goPath.toString()))
                                goInfo = jsonSlurper.parse(new File(goInfoPath.toString()))

                                def goIds = ipr2go[interproKey]
                                def goTerms = goIds.collect { goId ->
                                    goXref = [
                                        "name": goInfo[goId][0],
                                        "databaseName": "GO",
                                        "category": GO_PATTERN[goInfo[goId][1]],
                                        "id": goId
                                    ]
                                    GoXrefs goXrefObj = GoXrefs.fromMap(goXref)
                                    matchObject.signature.entry.addGoXrefs(goXrefObj)
                                }
                            } catch (java.lang.NullPointerException e) {
                                // pass
                            }
                        }

                        if (params.pathways) {
                            try {
                                String ipr2paPath = "${data_dir}/${params.xRefsConfig.pathways}.ipr.json"
                                String paInfoPath = "${data_dir}/${params.xRefsConfig.pathways}.json"
                                ipr2pa = jsonSlurper.parse(new File(ipr2paPath.toString()))
                                paInfo = jsonSlurper.parse(new File(paInfoPath.toString()))
                                def paIds = ipr2pa[interproKey]
                                def paTerms = paIds.collect { paId ->
                                    paXref = [
                                        "name": paInfo[paId][1],
                                        "databaseName": PA_PATTERN[paInfo[paId][0]],
                                        "id": paId
                                    ]
                                    PathwayXrefs paXrefObj = PathwayXrefs.fromMap(paXref)
                                    matchObject.signature.entry.addPathwayXrefs(paXrefObj)
                                }
                            } catch (java.lang.NullPointerException e) {
                                // pass
                            }
                        }
                    }
                }

                if (memberDB == "panther") {
                    String accSubfamily = matchObject.treegrafter.subfamilyAccession
                    matchObject.treegrafter.subfamilyAccession = accSubfamily
                    if (entriesInfo[accSubfamily]) {
                        matchObject.treegrafter.subfamilyName = entrieInfo[accSubfamily]["name"]
                    }
                }
                return [(matchId): matchObject]
            }]
        }
        def outputFilePath = task.workDir.resolve("matches2xrefs.json")
        def json = JsonOutput.toJson(matches)
        new File(outputFilePath.toString()).write(json)
    }
}
