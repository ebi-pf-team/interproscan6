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

    def (ipr2go, goInfo, ipr2pa, paInfo) = [null, null, null, null]
    if (params.goterms) {
        (ipr2go, goInfo) = loadXRefFiles("goterms", dataDir)
    }
    if (params.pathways) {
        (ipr2pa, paInfo) = loadXRefFiles("pathways", dataDir)
    }

    JsonSlurper jsonSlurper = new JsonSlurper()

    Map<String, Map<String, Match>> aggregatedMatches = [:]

    if ("${apps}".contains('panther')) {
        String paintAnnoDir = "${dataDir}/${params.appsConfig.paint}"
    }

    matchesEntries = membersMatches.each { matchesPath  ->
        memberDB = matchesPath.toString().split("/").last().split("\\.")[0]
        String memberRelease = null
        String library = null
        def libraryRelease = entries.databases.find { key, value ->
            InterProScan.standardiseMemberDB(key) == memberDB
        }
        try {
            memberRelease = libraryRelease.value
            library = libraryRelease.key
        } catch (java.lang.NullPointerException e) {
            if (memberDB == "cathfunfam") { // TODO: standardise FunFam on entries file and remove this block
                library = "cathfunfam"
            } else {
                throw new Exception("No library release found for ${memberDB} on entries data")
            }
        }
        SignatureLibraryRelease sigLibRelease = new SignatureLibraryRelease(
            library,
            memberRelease)

        def matches = jsonSlurper.parse(matchesPath).collectEntries { seqId, matches ->
            [(seqId): matches.collectEntries { rawModelAccession, match ->
                Match matchObject = Match.fromMap(match)
                def modelAccession = matchObject.modelAccession.split("\\.")[0]
                if (memberDB in ["cathgene3d", "cathfunfam", "superfamily"]) {
                    modelAccession = matchObject.signature.accession
                }

                if (!matchObject.signature) {
                    matchObject.signature = new Signature(modelAccession, sigLibRelease)
                } else {
                    matchObject.signature.signatureLibraryRelease = sigLibRelease
                }

                if (memberDB == "panther") {
                    String paintAnnPath = "${dataDir}/${params.appsConfig.panther.paint}/${matchObject.signature.accession}.json"
                    File paintAnnotationFile = new File(paintAnnPath)
                    // not every signature will have a paint annotation file match
                    if (paintAnnotationFile.exists()) {
                        def paintAnnotationsContent = jsonSlurper.parse(paintAnnotationFile)
                        String nodeId = matchObject.treegrafter.ancestralNodeID
                        def nodeData = paintAnnotationsContent[nodeId]
                        matchObject.treegrafter.subfamilyAccession = nodeData[0]
                        matchObject.treegrafter.proteinClass = nodeData[2]
                        matchObject.treegrafter.graftPoint = nodeData[3]
                    }
                }

                Map entry = entries['entries'][modelAccession] ?: entries['entries'][rawModelAccession]
                if (entry) {
                    matchObject.signature.name = entry["name"]
                    matchObject.signature.description = entry["description"]

                    if (entry['representative']) {
                        RepresentativeInfo representativeInfo = new RepresentativeInfo(
                            entry['representative']["type"],
                            entry['representative']["index"]
                        )
                        matchObject.representativeInfo = representativeInfo
                    }

                    def interproKey = entry['integrated']
                    def entryInfo = entries['entries'].get(interproKey)
                    if (entryInfo) {
                        Entry entryDataObj = new Entry(
                            interproKey,
                            entryInfo["name"],
                            entryInfo["description"],
                            entryInfo["type"]
                        )
                        matchObject.signature.entry = entryDataObj
                    }

                    if (params.goterms) {
                        try {
                            def goIds = ipr2go[interproKey]
                            def goTerms = goIds.collect { goId ->
                                GoXRefs goXRefObj = new GoXRefs(
                                    goInfo[goId][0],
                                    "GO",
                                    GO_PATTERN[goInfo[goId][1]],
                                    goId
                                )
                                matchObject.signature.entry.addGoXRefs(goXRefObj)
                            }
                        } catch (java.lang.NullPointerException e) {
                            // pass if no GO Terms found for the current entry
                        }
                    }
                    if (params.pathways) {
                        try {
                            def paIds = ipr2pa[interproKey]
                            def paTerms = paIds.collect { paId ->
                                PathwayXRefs paXRefObj = new PathwayXRefs(
                                    paInfo[paId][1],
                                    PA_PATTERN[paInfo[paId][0]],
                                    paId)
                                matchObject.signature.entry.addPathwayXRefs(paXRefObj)
                            }
                        } catch (java.lang.NullPointerException e) {
                            // pass if no Pathways found for the current entry
                        }
                    }
                }

                if (memberDB == "panther") {
                    accSubfamily = matchObject.signature.accession
                    if (entries['entries'][accSubfamily]) {
                        matchObject.treegrafter.subfamilyName = entries['entries'][accSubfamily]["name"]
                        matchObject.treegrafter.subfamilyDescription = entries['entries'][accSubfamily]["description"]
                    }
                }
                return [(rawModelAccession): matchObject]
            }]
        }
        matches.each { seqId, seqMatches ->
            if (aggregatedMatches.containsKey(seqId)) {
                seqMatches.each { rawModelAccession, matchObject ->
                    if (!aggregatedMatches[seqId].containsKey(rawModelAccession)) {
                        aggregatedMatches[seqId][rawModelAccession] = matchObject
                    }
                }
            } else {
                aggregatedMatches[seqId] = seqMatches
            }
        }
    }
    def outputFilePath = task.workDir.resolve("matches2xrefs.json")
    def json = JsonOutput.toJson(aggregatedMatches)
    new File(outputFilePath.toString()).write(json)
}

def loadXRefFiles(String xrefType, dataDir) {
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
