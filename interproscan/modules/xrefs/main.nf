import com.fasterxml.jackson.core.JsonFactory
import com.fasterxml.jackson.core.JsonParser
import com.fasterxml.jackson.core.JsonToken
import com.fasterxml.jackson.core.JsonGenerator
import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.core.JsonEncoding

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
    label 'small'

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
    def outputFilePath = task.workDir.resolve("matches2xrefs.json")
    String entriesPath = "${dataDir}/${entriesFile}"
    def (entries, ipr2go, goInfo, ipr2pa, paInfo) = [null, null, null, null, null]

    ObjectMapper mapper = new ObjectMapper()
    JsonFactory factory = new JsonFactory()
    JsonGenerator generator = factory.createGenerator(new File(outputFilePath.toString()), JsonEncoding.UTF8)
    generator.writeStartObject()

    if (!dataDir.toString().trim().isEmpty()) { // datadir doesn't need to be provided when only running members with no InterPro data
        File entriesJson = new File(entriesPath.toString())
        entries = mapper.readValue(entriesJson, Map)
        if (addGoterms) {
            (ipr2go, goInfo) = loadXRefFiles(gotermFilePrefix, dataDir)
        }
        if (addPathways) {
            (ipr2pa, paInfo) = loadXRefFiles(pathwaysFilePrefix, dataDir)
        }
    }

    matchesEntries = membersMatches.each { matchesPath  ->
        File matchFile = new File(matchesPath.toString())
        JsonParser parser = factory.createParser(matchFile)
        assert parser.nextToken() == JsonToken.START_OBJECT

        while (parser.nextToken() != JsonToken.END_OBJECT) {
            String seqId = parser.getCurrentName()
            parser.nextToken()
            generator.writeFieldName(seqId)
            generator.writeStartObject()

            while (parser.nextToken() != JsonToken.END_OBJECT) {
                String modelAccession = parser.getCurrentName()
                parser.nextToken()
                def match = mapper.readValue(parser, Map)
                Match matchObject = Match.fromMap(match)

                if (!entries) {
                    generator.writeFieldName(modelAccession)
                    generator.writeObject(matchObject)
                    continue
                }
                // null check needed for cases that signature still not created on match object (e.g. hmmer3 members)
                String entrySignatureKey = matchObject.signature?.accession ?: matchObject.modelAccession
                def signatureInfo = entries["entries"][entrySignatureKey] ?: entries["entries"][modelAccession]
                String memberDB = matchObject.signature?.signatureLibraryRelease?.library
                String memberRelease = null

                if (!memberDB) {
                    if (signatureInfo) {
                        memberDB = signatureInfo["database"]
                        memberRelease = entries["databases"][memberDB]
                    }
                    // PIRSR is the only memberDB that doesn't have entries associated on entries.json data file
                    else if (modelAccession.startsWith("PIRSR")) {
                        memberDB = "PIRSR"
                        memberRelease = entries["databases"][memberDB]
                    }
                    SignatureLibraryRelease sigLibRelease = new SignatureLibraryRelease(memberDB, memberRelease)
                    if (!matchObject.signature) {
                        matchObject.signature = new Signature(modelAccession, sigLibRelease)
                    } else if (!matchObject.signature.signatureLibraryRelease) {
                        matchObject.signature.signatureLibraryRelease = sigLibRelease
                    }
                }

                if (memberDB == "PANTHER" && matchObject.treegrafter.ancestralNodeID != null) {
                    String paintAnnPath = "${dataDir}/${paintAnnoDir}/${matchObject.signature.accession}.json"
                    File paintAnnotationFile = new File(paintAnnPath)
                    // not every signature will have a paint annotation file match
                    if (paintAnnotationFile.exists()) {
                        def paintAnnotationsContent = mapper.readValue(paintAnnotationFile, Map)
                        String nodeId = matchObject.treegrafter.ancestralNodeID
                        def nodeData = paintAnnotationsContent[nodeId]
                        matchObject.treegrafter.subfamilyAccession = nodeData[0]
                        matchObject.treegrafter.proteinClass = nodeData[2]
                        matchObject.treegrafter.graftPoint = nodeData[3]
                    }
                }

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
                    }

                    if (interproAcc && addGoterms) {
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
                    if (interproAcc && addPathways) {
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

                if (memberDB == "PANTHER") {
                    accSubfamily = matchObject.signature.accession
                    if (entries["entries"][accSubfamily]) {
                        matchObject.treegrafter.subfamilyName = entries["entries"][accSubfamily]["name"]
                        matchObject.treegrafter.subfamilyDescription = entries["entries"][accSubfamily]["description"]
                    }
                }

                generator.writeFieldName(modelAccession)
                mapper.writeValue(generator, matchObject)            
            }

            generator.writeEndObject()
        }

        parser.close()
    }

    generator.writeEndObject()
    generator.close()
}

def loadXRefFiles(xrefDir, dataDir) {
    String iprFilePath = "${dataDir}/${xrefDir}.ipr.json"
    String infoFilePath = "${dataDir}/${xrefDir}.json"
    File iprFile = new File(iprFilePath.toString())
    File infoFile = new File(infoFilePath.toString())
    ObjectMapper mapper = new ObjectMapper()

    if (!iprFile.exists()) { throw new FileNotFoundException("${iprFilePath} file not found") }
    if (!infoFile.exists()) { throw new FileNotFoundException("${infoFile} file not found") }

    try {
        def iprData = mapper.readValue(iprFile, Map)
        def infoData = mapper.readValue(infoFile, Map)
        return [iprData, infoData]
    } catch (Exception e) {
        throw new Exception("Error parsing goterms/pathways files: ${e}")
    }
}

