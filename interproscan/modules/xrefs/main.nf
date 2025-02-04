import com.fasterxml.jackson.core.JsonGenerator
import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.databind.JsonNode
import com.fasterxml.jackson.databind.SerializationFeature

// Define mapping patterns
def GO_PATTERN = ["P": "BIOLOGICAL_PROCESS", "C": "CELLULAR_COMPONENT", "F": "MOLECULAR_FUNCTION"]
def PA_PATTERN = ["t": "MetaCyc", "w": "UniPathway", "k": "KEGG", "r": "Reactome"]

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
        entries = JsonReader.load(entriesPath, jacksonMapper)
        if (addGoterms) (ipr2go, goInfo) = loadXRefFiles(gotermFilePrefix, dataDir, jacksonMapper)
        if (addPathways) (ipr2pa, paInfo) = loadXRefFiles(pathwaysFilePrefix, dataDir, jacksonMapper)
    }

    // Output JSON writing
    String outputFilePath = task.workDir.resolve("matches2xrefs.json")
    JsonWriter.stream(outputFilePath, jacksonMapper) { JsonGenerator jsonGenerator ->
        jsonGenerator.writeStartObject()
        membersMatches.each { matchesPath ->
            File file = new File(matchesPath.toString())
            JsonReader.streamJson(matchesPath.toString(), jacksonMapper) { String seqId, JsonNode matches ->
                jsonGenerator.writeFieldName(seqId)
                jsonGenerator.writeStartObject()
                
                matches.fields().each { Map.Entry<String, JsonNode> entry ->
                    String modelAcc = entry.key   // Extract the modelAcc (key)
                    Match match = Match.fromJsonNode(entry.value)
                    
                        if (!entries) {
                            jsonGenerator.writeFieldName(modelAcc)
                            JsonWriter.writeMap(jsonGenerator, match)
                        } else {
                            String signatureAcc = match.signature.accession
                            def signatureInfo = entries["entries"][signatureAcc] ?: entries["entries"][modelAcc]

                            // Update library version
                            if (!match.signature.signatureLibraryRelease.version && signatureInfo) {
                                match.signature.signatureLibraryRelease.version = entries["databases"][signatureInfo["database"]]
                            }

                            // Handle PANTHER data
                            if (match.signature.signatureLibraryRelease.library == "PANTHER") {
                                updatePantherData(match, dataDir, paintAnnoDir, signatureAcc, entries)
                            }

                            // Update signature info
                            if (signatureInfo) {
                                match.signature.name = signatureInfo["name"]
                                match.signature.description = signatureInfo["description"]
                                if (signatureInfo["representative"]) {
                                    match.representativeInfo = new RepresentativeInfo(
                                        signatureInfo["representative"]["type"],
                                        signatureInfo["representative"]["index"]
                                    )
                                }

                                // Handle InterPro data
                                def interproAcc = signatureInfo["integrated"]
                                if (interproAcc) {
                                    def entryInfo = entries["entries"].get(interproAcc)
                                    assert entryInfo != null
                                    match.signature.entry = new Entry(
                                        interproAcc, entryInfo["name"], entryInfo["description"], entryInfo["type"]
                                    )
                                    addXRefs(match, interproAcc, ipr2go, goInfo, ipr2pa, paInfo)
                                }
                            }
                            // Write out the model Acc and updated Match object
                            jsonGenerator.writeFieldName(modelAcc)
                            JsonWriter.writeMap(jsonGenerator, match)
                        }  // end of if/else
                    jsonGenerator.writeEndObject()
                } // end of matches
            }  // end of Json reader / seq Id
        } // end of members matches
        jsonGenerator.writeEndObject()
    }  // end of Json writer stream
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
    if (ipr2go && goInfo && ipr2go[interproAcc]) {
        ipr2go[interproAcc].each { goId ->
            match.signature.entry.addGoXRefs(new GoXRefs(goInfo[goId][0], "GO", GO_PATTERN[goInfo[goId][1]], goId))
        }
    }
    if (ipr2pa && paInfo && ipr2pa[interproAcc]) {
        ipr2pa[interproAcc].each { paId ->
            match.signature.entry.addPathwayXRefs(new PathwayXRefs(paInfo[paId][1], PA_PATTERN[paInfo[paId][0]], paId))
        }
    }
}
