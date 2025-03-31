import com.fasterxml.jackson.core.JsonGenerator
import com.fasterxml.jackson.databind.ObjectMapper
import groovy.json.JsonOutput

process XREFS {
    label 'local', 'tiny'

    input:
    tuple val(meta), val(membersMatches)
    val apps
    val dataDir
    val xrefDir
    val memberDbReleases
    val entriesFile
    val gotermFilePrefix
    val pathwaysFilePrefix
    val addGoterms
    val addPathways
    val paintAnnoDir

    output:
    tuple val(meta), path("matches2xrefs.json")

    exec:
    // Load entries data. NOTE: The entries, go terms and pathway files are too large for JsonSlurper to handle.
    String entriesPath = "${xrefDir}/${entriesFile}"
    def (entries, ipr2go, goInfo, ipr2pa, paInfo) = [null, null, null, null, null]
    if (dataDir.toString().trim()) {  // datadir is not needed when exclusively running members with no interpro data
        File entriesJson = new File(entriesPath)
        entries = new ObjectMapper().readValue(entriesJson, Map)
        if (addGoterms) (ipr2go, goInfo) = loadXRefFiles(gotermFilePrefix, xrefDir)
        if (addPathways) (ipr2pa, paInfo) = loadXRefFiles(pathwaysFilePrefix, xrefDir)
    }

    def aggregatedMatches = [:]  // seqMd5: {modelAcc: Match} -- otherwise a seqMd5 will appear multiple times in the output
    membersMatches.each { matchesPath ->
        def matchesFileMap = new ObjectMapper().readValue(new File(matchesPath.toString()), Map)
        matchesFileMap.each { String seqMd5, Map matches ->
            def seqEntry = aggregatedMatches.computeIfAbsent(seqMd5, { [:] } )

            matches.each { modelAcc, match ->
                match = Match.fromMap(match)  // convert Map to Match object

                if (!entries || entries == null) {  // no data to update, update match in aggregatedMatches
                    seqEntry[modelAcc] = match
                } else {
                    String signatureAcc = match.signature.accession
                    def signatureInfo = entries[signatureAcc] ?: entries[modelAcc]

                    // Update library version
                    def version = (match.signature.signatureLibraryRelease.version == "null") ? null : match.signature.signatureLibraryRelease.version
                    if (!version && signatureInfo != null) {
                        def memberName = InterPro.formatMemberDbName(signatureInfo["database"])
                        match.signature.signatureLibraryRelease.version = memberDbReleases[memberName]
                    }

                    // Handle PANTHER data
                    if (match.signature.signatureLibraryRelease.library == "PANTHER") {
                        updatePantherData(match, paintAnnoDir, signatureAcc, entries)
                    }

                    // Update signature info
                    if (signatureInfo != null) {
                        match.signature.name = signatureInfo["name"]
                        match.signature.description = signatureInfo["description"]
                        String sigType = signatureInfo["type"]
                        match.signature.setType(sigType)

                        if (signatureInfo["representative"] != null) {
                            match.representativeInfo = new RepresentativeInfo(
                                signatureInfo["representative"]["type"],
                                signatureInfo["representative"]["index"]
                            )
                        }

                        // Handle InterPro data
                        String interproAcc = signatureInfo["integrated"]
                        if (interproAcc != null) {
                             def entryInfo = entries[interproAcc]
                             assert entryInfo != null
                             match.signature.entry = new Entry(
                                 interproAcc, entryInfo["name"], entryInfo["description"], entryInfo["type"]
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

def loadXRefFiles(xrefFileName, xrefDir) {
    def iprFilePath = new File("${xrefDir}/${xrefFileName}.ipr.json")
    def infoFilePath = new File("${xrefDir}/${xrefFileName}.json")
    return [
        new ObjectMapper().readValue(iprFilePath, Map),
        new ObjectMapper().readValue(infoFilePath, Map)
    ]
}

def updatePantherData(def match, def paintAnnoDir, def signatureAcc, def entries) {
    File paintAnnotationFile = new File("${paintAnnoDir}/${signatureAcc}.json")
    if (paintAnnotationFile.exists()) {
        def paintAnnotationsContent = new ObjectMapper().readValue(paintAnnotationFile, Map)
        String nodeId = match.treegrafter.ancestralNodeID
        def nodeData = paintAnnotationsContent[nodeId]
        if (nodeData != null) {
            match.treegrafter.subfamilyAccession = (nodeData[0] == "null") ? null : nodeData[0]
            match.treegrafter.proteinClass = (nodeData[2] == "null") ? null : nodeData[2]
            match.treegrafter.graftPoint = (nodeData[3] == "null") ? null : nodeData[3]
        }
    }
    if (entries[signatureAcc]) {
        match.treegrafter.subfamilyName = entries[signatureAcc]["name"]
        match.treegrafter.subfamilyDescription = entries[signatureAcc]["description"]
    }
}

def addXRefs(Match match, String interproAcc, def ipr2go, def goInfo, def ipr2pa, def paInfo) {
    Map<String,String> GO_PATTERN = ["P": "BIOLOGICAL_PROCESS", "C": "CELLULAR_COMPONENT", "F": "MOLECULAR_FUNCTION"]
    Map<String,String> PA_PATTERN = ["t": "MetaCyc", "w": "UniPathway", "k": "KEGG", "r": "Reactome"]
    if (ipr2go != null && goInfo != null && ipr2go[interproAcc] != null) {
        ipr2go[interproAcc].each { goId ->
            goId = goId
            match.signature.entry.addGoXRefs(
                new GoXRefs(goInfo["terms"][goId][0], "GO", GO_PATTERN[goInfo["terms"][goId][1]], goId)
            )
        }
    }
    if (ipr2pa != null && paInfo != null && ipr2pa[interproAcc] != null) {
        ipr2pa[interproAcc].each { paId ->
            paId = paId
            match.signature.entry.addPathwayXRefs(
                new PathwayXRefs(paInfo[paId][1], PA_PATTERN[paInfo[paId][0]], paId)
            )
        }
    }
}
