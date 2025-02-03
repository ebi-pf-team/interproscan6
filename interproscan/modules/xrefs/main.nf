import com.fasterxml.jackson.core.JsonGenerator
import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.databind.JsonNode
import com.fasterxml.jackson.databind.SerializationFeature

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
    // Build a single mapper for all readers and writers to save memory
    ObjectMapper jacksonMapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT)

    // Load entries data
    String entriesPath = "${dataDir}/${entriesFile}"
    def (entries, ipr2go, goInfo, ipr2pa, paInfo) = [null, null, null, null, null]
    if (!dataDir.toString().trim().isEmpty()) { // datadir doesn't need to be provided when only running members with no InterPro data
        File entriesJson = new File(entriesPath.toString())
        entries = JsonReader.load(entriesPath.toString(), jacksonMapper)
        if (addGoterms) {
            (ipr2go, goInfo) = loadXRefFiles(gotermFilePrefix, dataDir, jacksonMapper)
        }
        if (addPathways) {
            (ipr2pa, paInfo) = loadXRefFiles(pathwaysFilePrefix, dataDir, jacksonMapper)
        }
    }

    // Stream writing the JSON file of matches with added XREFs
    String outputFilePath = task.workDir.resolve("matches2xrefs.json")
    JsonWriter.stream(outputFilePath.toString(), jacksonMapper) { JsonGenerator jsonGenerator ->
        jsonGenerator.writeStartObject()
        membersMatches.each { matchesPath  ->
            // matchesPath is keyed by the prot seq id and valued by a map [modelAcc: [match data = match object]]
            JsonReader.stream(matchesPath.toString(), jacksonMapper) { JsonNode seq ->
                seq.fields().each { entry ->
                    String seqId = entry.key
                    jsonGenerator.writeFieldName(seqId)
                    jsonGenerator.writeStartObject()  // start new Json object, e.g. nested Map

                    JsonNode matches = entry.value  // matches for this seq: [modelAccession: [match data]]
                    matches.fields().each { match
                        String modelAcc = match.key
                        Match matchObject = Match.fromJsonNode(match.value)

                        if (!entries) {
                            jsonGenerator.writeFieldName(modelAcc)
                            JsonWriter.writeMap(jsonGenerator, matchObject)
                        } else {
                            // Retrieve the relevant signature info from the match and entries
                            String signatureAcc = matchObject.signature.accession
                            def signatureInfo = entries["entries"][signatureAcc] ?: entries["entries"][modelAcc]

                            // Update the library version if the sig/model is found in the JSON
                            String memberRelease = matchObject.signature.signatureLibraryRelease.version
                            if (!memberRelease && signatureInfo) {
                                matchObject.signature.signatureLibraryRelease.version = entries["databases"][signatureInfo["database"]]
                            }

                            // PANTHER: Add PAINT annotation data and subfamily name and description
                            if (matchObject.signature.signatureLibraryRelease.library == "PANTHER") {
                                if (matchObject.treegrafter.ancestralNodeID != null) {
                                    String paintAnnPath = "${dataDir}/${paintAnnoDir}/${signatureAcc}.json"
                                    File paintAnnotationFile = new File(paintAnnPath)
                                    // not every signature will have a PAINT annotation file match
                                    if (paintAnnotationFile.exists()) {
                                        def paintAnnotationsContent = JsonProcessor.jsonToMap(paintAnnotationFile)
                                        String nodeId = matchObject.treegrafter.ancestralNodeID
                                        def nodeData = paintAnnotationsContent[nodeId]
                                        matchObject.treegrafter.subfamilyAccession = nodeData[0]
                                        matchObject.treegrafter.proteinClass = nodeData[2]
                                        matchObject.treegrafter.graftPoint = nodeData[3]
                                    }
                                }

                                if (entries["entries"][signatureAcc]) {
                                    matchObject.treegrafter.subfamilyName = entries["entries"][signatureAcc]["name"]
                                    matchObject.treegrafter.subfamilyDescription = entries["entries"][signatureAcc]["description"]
                                }
                            }

                            // Update signature info if there is any
                            if (signatureInfo) {
                                matchObject.signature.name = signatureInfo["name"]
                                matchObject.signature.description = signatureInfo["description"]

                                if (signatureInfo["representative"]) {
                                    RepresentativeInfo representativeInfo = new RepresentativeInfo(
                                        signatureInfo["representative"]["type"],
                                        signatureInfo["representative"]["index"]
                                    )
                                    matchObject.representativeInfo = representativeInfo
                                }

                                // add InterPro data if there match is intergrated into any InterPro entries
                                def interproAcc = signatureInfo["integrated"]
                                if (interproAcc) {
                                    def entryInfo = entries["entries"].get(interproAcc)
                                    assert entryInfo != null
                                    Entry entryDataObj = new Entry(
                                        interproAcc,
                                        entryInfo["name"],
                                        entryInfo["description"],
                                        entryInfo["type"]
                                    )
                                    matchObject.signature.entry = entryDataObj

                                    if (addGoterms) {
                                        def goIds = ipr2go[interproAcc]
                                        if (goIds) {
                                            def goTerms = goIds.collect { goId ->
                                                GoXRefs goXRefObj = new GoXRefs(
                                                    goInfo[goId][0],
                                                    "GO",
                                                    GO_PATTERN[goInfo[goId][1]],
                                                    goId
                                                )
                                                matchObject.signature.entry.addGoXRefs(goXRefObj)
                                            }
                                        }
                                    }

                                    if (addPathways) {
                                        def paIds = ipr2pa[interproAcc]
                                        if (paIds) {
                                            def paTerms = paIds.collect { paId ->
                                                PathwayXRefs paXRefObj = new PathwayXRefs(
                                                    paInfo[paId][1],
                                                    PA_PATTERN[paInfo[paId][0]],
                                                    paId)
                                                matchObject.signature.entry.addPathwayXRefs(paXRefObj)
                                            }
                                        }
                                    }
                                }
                            }

                            // Write out the model Acc and updated Match object
                            jsonGenerator.writeFieldName(modelAcc)
                            JsonWriter.writeMap(jsonGenerator, matchObject)
                        }
                    }
                    jsonGenerator.writeEndObject()  // end of seqId
                }
            }
        }
        jsonGenerator.writeEndObject()  // end of whole object
    }
}

def loadXRefFiles(xrefDir, dataDir, jacksonMapper) {
    String iprFilePath = "${dataDir}/${xrefDir}.ipr.json"
    String infoFilePath = "${dataDir}/${xrefDir}.json"
    if (!new File(iprFilePath.toString()).exists()) { throw new FileNotFoundException("${iprFilePath} file not found") }
    if (!new File(infoFilePath.toString()).exists()) { throw new FileNotFoundException("${infoFilePath} file not found") }

    try {
        def iprData = JsonReader.load(iprFilePath, jacksonMapper)
        def infoData = JsonReader.load(infoFilePath, jacksonMapper)
        return [iprData, infoData]
    } catch (Exception e) {
        throw new Exception("Error parsing goterms/pathways files: ${e}")
    }
}
