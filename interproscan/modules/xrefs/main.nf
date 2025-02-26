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
        println "entries class = ${entries.getClass()}"
        if (addGoterms) (ipr2go, goInfo) = loadXRefFiles(gotermFilePrefix, dataDir, jacksonMapper)
        if (addPathways) (ipr2pa, paInfo) = loadXRefFiles(pathwaysFilePrefix, dataDir, jacksonMapper)
    }

    def aggregatedMatches = [:]  // seqId: {modelAcc: Match} -- otherwise a seqId will appear multiple times in the output
    membersMatches.each { matchesPath ->
        File file = new File(matchesPath.toString())
        JsonReader.streamJson(matchesPath.toString(), jacksonMapper) { String seqId, JsonNode matches ->
            assert seqId != null : "Error: seqId is null in $matchesPath"
            seqEntry = aggregatedMatches.computeIfAbsent(seqId, { [:] } )

            matches.fields().each { Map.Entry<String, JsonNode> entry ->
                String modelAcc = entry.key   // Extract the modelAcc (key)
                Match match = Match.fromJsonNode(entry.value)

                if (entries != null) {  // no data to update, update match in aggregatedMatches
                    seqEntry[modelAcc] = match
                } else {
                    String signatureAcc = match.signature.accession
                    def signatureInfo = entries["entries"][signatureAcc] ?: entries["entries"][modelAcc]

                    // Update library version
                    def version = (match.signature.signatureLibraryRelease.version == "null") ? null : match.signature.signatureLibraryRelease.version
                    if (!version && signatureInfo != null) {
                        match.signature.signatureLibraryRelease.version = entries["databases"][signatureInfo["database"].asText()].asText()
                    }

                    // Handle PANTHER data
                    if (match.signature.signatureLibraryRelease.library == "PANTHER") {
                        updatePantherData(match, dataDir, paintAnnoDir, signatureAcc, entries)
                    }

                    // Update signature info
                    if (signatureInfo != null) {
                        match.signature.name = signatureInfo["name"].asText().replace('\"',"")
                        match.signature.description = signatureInfo["description"].asText().replace('\"',"")
                        String sigType = signatureInfo["type"].asText().replace('\"',"")
                        matchObject.signature.setType(sigType)
                        if (signatureInfo["representative"] != null) { // if(var) does not work on JsonNodes, need explicit falsey check
                            match.representativeInfo = new RepresentativeInfo(
                                signatureInfo["representative"].get("type").asText(),
                                signatureInfo["representative"].get("index").intValue()
                            )
                        }

                        // Handle InterPro data
                        String interproAcc = signatureInfo["integrated"].asText()  // defaults to a TextNode, convert to String
                        if (interproAcc != "null") {
                             def entryInfo = entries["entries"][interproAcc]
                             assert entryInfo != null
                             match.signature.entry = new Entry(
                                 interproAcc, entryInfo.get("name").asText(), entryInfo.get("description").asText(), entryInfo.get("type").asText()
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

    // Stream writing the output JSON file
    String outputFilePath = task.workDir.resolve("matches2xrefs.json")
    JsonWriter.writeMaptoFile(outputFilePath.toString(), jacksonMapper, aggregatedMatches)
}

def loadXRefFiles(xrefDir, dataDir, jacksonMapper) {
    String iprFilePath = "${dataDir}/${xrefDir}.ipr.json"
    String infoFilePath = "${dataDir}/${xrefDir}.json"
    if (!new File(iprFilePath).exists() || !new File(infoFilePath).exists()) {
        throw new FileNotFoundException("XRef files missing: ${iprFilePath}, ${infoFilePath}")
    }
    return [JsonReader.load(iprFilePath, jacksonMapper), JsonReader.load(infoFilePath, jacksonMapper)]
}

def updatePantherData(Match match, String dataDir, String paintAnnoDir, String signatureAcc, Map entries) {
    String paintAnnPath = "${dataDir}/${paintAnnoDir}/${signatureAcc}.json"
    File paintAnnotationFile = new File(paintAnnPath)
    if (paintAnnotationFile.exists()) {
        def paintAnnotationsContent = JsonProcessor.jsonToMap(paintAnnotationFile)
        String nodeId = match.treegrafter.ancestralNodeID
        def nodeData = paintAnnotationsContent[nodeId]
        if (nodeData) {
            match.treegrafter.subfamilyAccession = nodeData[0]
            match.treegrafter.proteinClass = nodeData[2]
            match.treegrafter.graftPoint = nodeData[3]
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
            goId = goId.asText().replaceAll('"', '')
            match.signature.entry.addGoXRefs(
                new GoXRefs(goInfo.get("terms").get(goId).get(0).asText(), "GO", GO_PATTERN[goInfo.get("terms").get(goId).get(1).asText()], goId)
            )
        }
    }
    if (ipr2pa != null && paInfo != null && ipr2pa[interproAcc] != null) {
        ipr2pa[interproAcc].each { paId ->
            paId = paId.asText().replaceAll('"', '')
            match.signature.entry.addPathwayXRefs(new PathwayXRefs(
                paInfo.get(paId).get(1).asText(), PA_PATTERN[paInfo.get(paId).get(0).asText()], paId
            ))
        }
    }
}