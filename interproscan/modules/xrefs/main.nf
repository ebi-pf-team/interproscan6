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
    val dataDir

    output:
    tuple val(meta), path("matches2xrefs.json")

    exec:
    String entriesPath = "${dataDir}/${params.xRefsConfig.entries}"
    File entriesJson = new File(entriesPath.toString())
    def entries = new ObjectMapper().readValue(entriesJson, Map)

    def (ipr2go, goInfo) = loadXrefFiles("goterms", dataDir)
    def (ipr2pa, paInfo) = loadXrefFiles("pathways", dataDir)

    def jsonOutput = []
    JsonSlurper jsonSlurper = new JsonSlurper()
    matchesEntries = membersMatches.each { matchesPath  ->
        memberDB = matchesPath.toString().split("/").last().split("\\.")[0]
        SignatureLibraryRelease sigLibRelease = new SignatureLibraryRelease(memberDB, entries['databases'][memberDB])
        def matches = jsonSlurper.parse(matchesPath).collectEntries { seqId, jsonMatches ->
            [(seqId): jsonMatches.collectEntries { matchId, jsonMatch ->
                Match matchObject = Match.fromMap(jsonMatch)
                def entriesInfo = entries['entries']
                String accId = matchObject.modelAccession.split("\\.")[0]
                Signature signatureObject = new Signature(accId, "", "", sigLibRelease, null)
                matchObject.signature = signatureObject

                if (memberDB == "panther") {
                    String sigAcc = matchObject.signature.accession
                    String paintAnnPath = "${params.appsConfig.panther.paintAnnDir}/${sigAcc}.json"
                    File paintAnnotationFile = new File(paintAnnPath.toString())
                    if (paintAnnotationFile.exists()) {
                        def paintAnnotationsContent = jsonSlurper.parse(paintAnnotationFile)
                        String nodeId = matchObject.treegrafter.ancestralNodeID
                        def nodeData = paintAnnotationsContent[nodeId]
                        matchObject.treegrafter.proteinClass = nodeData[2]
                        matchObject.treegrafter.graftPoint = nodeData[3]
                    }
                }

                def entry = entriesInfo[accId] ?: entriesInfo[matchId]
                if (entry) {
                    matchObject.signature.name = entry["name"]
                    matchObject.signature.description = entry["description"]

                    def representativeFlag = null
                    if (entry['representative']) {
                        RepresentativeInfo representativeInfo = new RepresentativeInfo(
                            entry['representative']["type"],
                            rank = entry['representative']["index"]
                        )
                        matchObject.representativeInfo = representativeInfo
                    }

                    def interproKey = entry['integrated']
                    def entryInfo = entriesInfo.get(interproKey)
                    if (entryInfo) {
                        Entry entryDataObj = new Entry(
                            interproKey,
                            entryInfo["name"],
                            entryInfo["description"],
                            entryInfo["type"]
                        )
                        matchObject.signature.entry = entryDataObj
                    }

                    if (params.goterms || params.pathways) {
                        interproKey = entry["integrated"]
                        if (params.goterms) {
                            try {
                                def goIds = ipr2go[interproKey]
                                def goTerms = goIds.collect { goId ->
                                    GoXrefs goXrefObj = new GoXrefs(
                                        goInfo[goId][0],
                                        "GO",
                                        GO_PATTERN[goInfo[goId][1]],
                                        goId
                                    )
                                    matchObject.signature.entry.addGoXrefs(goXrefObj)
                                }
                            } catch (java.lang.NullPointerException e) {
                                // pass
                            }
                        }

                        if (params.pathways) {
                            try {
                                def paIds = ipr2pa[interproKey]
                                def paTerms = paIds.collect { paId ->
                                    PathwayXrefs paXrefObj = new PathwayXrefs(
                                        paInfo[paId][1],
                                        PA_PATTERN[paInfo[paId][0]],
                                        paId)
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
        jsonOutput.add(matches)
    }
    def outputFilePath = task.workDir.resolve("matches2xrefs.json")
    def json = JsonOutput.toJson(jsonOutput)
    new File(outputFilePath.toString()).write(json)
}

def loadXrefFiles(String xrefType, dataDir) {
    JsonSlurper jsonSlurper = new JsonSlurper()

    String iprFilePath = "${dataDir}/${params.xRefsConfig[xrefType]}.ipr.json"
    String infoFilePath = "${dataDir}/${params.xRefsConfig[xrefType]}.json"
    File iprFile = new File(iprFilePath.toString())
    File infoFile = new File(infoFilePath.toString())

    if (!iprFile.exists()) { throw new FileNotFoundException("${iprFilePath} file not found") }
    if (!infoFile.exists()) { throw new FileNotFoundException("${infoFile} file not found") }

    try {
        def iprData = jsonSlurper.parse(iprFile)
        def infoData = jsonSlurper.parse(infoFile)
        return [iprData, infoData]
    } catch (Exception e) {
        throw new Exception("Error parsing ${xrefType} files: ${e}")
    }
}
