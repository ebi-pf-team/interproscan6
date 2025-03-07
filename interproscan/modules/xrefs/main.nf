import com.fasterxml.jackson.core.JsonGenerator
import com.fasterxml.jackson.databind.JsonNode
import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.databind.SerializationFeature

process XREFS {
    label 'local'

    input:
    tuple val(meta), val(membersMatches)
    val apps
    val dataDir
    val entriesFile
    val gotermFilePrefix
    val pathwaysFilePrefix
    val addGoterms
    val addPathways
    val paintAnnoDir

    output:
    tuple val(meta), path("matches2xrefs.json")

    exec:
    ObjectMapper jacksonMapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT)

    // Load entries data
    String entriesPath = "${dataDir}/${entriesFile}"
    def (entries, ipr2go, goInfo, ipr2pa, paInfo) = [null, null, null, null, null]
    if (dataDir.toString().trim()) {  // datadir is not needed when exclusively running members with no interpro data
        File entriesJson = new File(entriesPath)
        entries = JsonReader.load(entriesPath, jacksonMapper) // returns
        if (addGoterms) (ipr2go, goInfo) = loadXRefFiles(gotermFilePrefix, dataDir, jacksonMapper)
        if (addPathways) (ipr2pa, paInfo) = loadXRefFiles(pathwaysFilePrefix, dataDir, jacksonMapper)
    }

    def aggregatedMatches = [:]  // seqMd5: {modelAcc: Match} -- otherwise a seqMd5 will appear multiple times in the output
    membersMatches.each { matchesPath ->
        File file = new File(matchesPath.toString())
        JsonReader.streamJson(matchesPath.toString(), jacksonMapper) { String seqMd5, JsonNode matches ->
            seqEntry = aggregatedMatches.computeIfAbsent(seqMd5, { [:] } )

            matches.fields().each { Map.Entry<String, JsonNode> entry ->
                String modelAcc = entry.key   // Extract the modelAcc (key)
                Match match = Match.fromJsonNode(entry.value)

                if (!entries || entries == null) {  // no data to update, update match in aggregatedMatches
                    seqEntry[modelAcc] = match
                } else {
                    String signatureAcc = match.signature.accession
                    def signatureInfo = entries["entries"][signatureAcc] ?: entries["entries"][modelAcc]

                    // Update library version
                    def version = (match.signature.signatureLibraryRelease.version == "null") ? null : match.signature.signatureLibraryRelease.version
                    if (!version && signatureInfo != null) {
                        match.signature.signatureLibraryRelease.version = JsonReader.asString(entries["databases"][JsonReader.asString(signatureInfo["database"])])
                    }

                    // Handle PANTHER data
                    if (match.signature.signatureLibraryRelease.library == "PANTHER") {
                        updatePantherData(match, dataDir, paintAnnoDir, signatureAcc, entries, jacksonMapper)
                    }

                    // Update signature info
                    if (signatureInfo != null) {  // signatureInfo is an objectNode !var does not work
                        match.signature.name = JsonReader.asString(signatureInfo["name"])
                        match.signature.description = JsonReader.asString(signatureInfo["description"])
                        String sigType = JsonReader.asString(signatureInfo["type"])
                        match.signature.setType(sigType)

                        if (signatureInfo["representative"] != null) { // if(var) does not work on JsonNodes, need explicit falsey check
                            match.representativeInfo = new RepresentativeInfo(
                                JsonReader.asString(signatureInfo["representative"].get("type"))
                                signatureInfo["representative"].get("index").intValue()
                            )
                        }

                        // Handle InterPro data
                        String interproAcc = JsonReader.asString(signatureInfo["integrated"])  // defaults to a TextNode, convert to String
                        if (interproAcc != null) {
                             def entryInfo = entries["entries"][interproAcc]
                             assert entryInfo != null
                             match.signature.entry = new Entry(
                                 interproAcc, JsonReader.asString(entryInfo.get("name")), JsonReader.asString(entryInfo.get("description")), JsonReader.asString(entryInfo.get("type"))
                             )
                             addXRefs(match, interproAcc, ipr2go, goInfo, ipr2pa, paInfo)
                        }
                    }
                    // Update the in aggregatedMatches Match object
                    seqEntry[modelAcc] = match
                }  // end of if/else
            } // end of matches
        }  // end of Json reader / seq Id
    } // end of members matches

    String outputFilePath = task.workDir.resolve("matches2xrefs.json")
    def json = JsonOutput.toJson(aggregatedMatches)
    new File(outputFilePath.toString()).write(json)
}

def loadXRefFiles(xrefDir, dataDir, jacksonMapper) {
    String iprFilePath = "${dataDir}/${xrefDir}.ipr.json"
    String infoFilePath = "${dataDir}/${xrefDir}.json"
    if (!new File(iprFilePath).exists() || !new File(infoFilePath).exists()) {
        throw new FileNotFoundException("XRef files missing: ${iprFilePath}, ${infoFilePath}")
    }
    return [JsonReader.load(iprFilePath, jacksonMapper), JsonReader.load(infoFilePath, jacksonMapper)]
}

def updatePantherData(def match, def dataDir, def paintAnnoDir, def signatureAcc, def entries, def jacksonMapper) {
    String paintAnnPath = "${dataDir}/${paintAnnoDir}/${signatureAcc}.json"
    File paintAnnotationFile = new File(paintAnnPath)
    if (paintAnnotationFile.exists()) {
        def paintAnnotationsContent = JsonReader.load(paintAnnotationFile.toString(), jacksonMapper)
        String nodeId = match.treegrafter.ancestralNodeID
        def nodeData = paintAnnotationsContent.get(nodeId)
        if (nodeData != null) {
            match.treegrafter.subfamilyAccession = (nodeData[0] == "null") ? null : nodeData[0]
            match.treegrafter.proteinClass = (nodeData[2] == "null") ? null : nodeData[2]
            match.treegrafter.graftPoint = (nodeData[3] == "null") ? null : nodeData[3]
        }
    }
    if (entries["entries"][signatureAcc]) {
        match.treegrafter.subfamilyName = entries["entries"][signatureAcc]["name"]
        match.treegrafter.subfamilyDescription = entries["entries"][signatureAcc]["description"]
    }
}

def addXRefs(Match match, String interproAcc, def ipr2go, def goInfo, def ipr2pa, def paInfo) {
    Map<String,String> GO_PATTERN = ["P": "BIOLOGICAL_PROCESS", "C": "CELLULAR_COMPONENT", "F": "MOLECULAR_FUNCTION"]
    Map<String,String> PA_PATTERN = ["t": "MetaCyc", "w": "UniPathway", "k": "KEGG", "r": "Reactome"]
    if (ipr2go != null && goInfo != null && ipr2go[interproAcc] != null) {
        ipr2go[interproAcc].each { goId ->
            goId = JsonReader.asString(goId)
            match.signature.entry.addGoXRefs(
                new GoXRefs(JsonReader.asString(goInfo.get("terms").get(goId).get(0)), "GO", GO_PATTERN[JsonReader.asString(goInfo.get("terms").get(goId).get(1))], goId)
            )
        }
    }
    if (ipr2pa != null && paInfo != null && ipr2pa[interproAcc] != null) {
        ipr2pa[interproAcc].each { paId ->
            paId = JsonReader.asString(paId)
            match.signature.entry.addPathwayXRefs(
                new PathwayXRefs(JsonReader.asString(paInfo.get(paId).get(1)), PA_PATTERN[JsonReader.asString(paInfo.get(paId).get(0))], paId)
            )
        }
    }
}