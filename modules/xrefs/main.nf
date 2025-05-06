import com.fasterxml.jackson.core.JsonGenerator
import com.fasterxml.jackson.databind.ObjectMapper
import groovy.json.JsonOutput

process XREFS {
    executor 'local'

    input:
    tuple val(meta), val(members_matches)
    val applications
    val db_releases
    val add_goterms
    val add_pathways
    val panther_paint_dir

    output:
    tuple val(meta), path("matches2xrefs.json")

    exec:
    def (databaseInfo, entries, ipr2go, goInfo, ipr2pa, paInfo) = [null, null, null, null, null, null]
    if (db_releases.interpro) {
        String interproDir = db_releases.interpro.dirpath.toString()
        String databasesPath = "${interproDir}/databases.json"
        File databasesJson = new File(databasesPath)
        databaseInfo = new ObjectMapper().readValue(databasesJson, Map)
        String entriesPath = "${interproDir}/entries.json"
        File entriesJson = new File(entriesPath)
        entries = new ObjectMapper().readValue(entriesJson, Map)
        if (add_goterms) (ipr2go, goInfo) = loadXRefFiles("${interproDir}/goterms")
        if (add_pathways) (ipr2pa, paInfo) = loadXRefFiles("${interproDir}/pathways")
    }

    def aggregatedMatches = [:]  // seqMd5: {modelAcc: Match} -- otherwise a seqMd5 will appear multiple times in the output
    members_matches.each { matchesPath ->
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
                        match.signature.signatureLibraryRelease.version = databaseInfo[signatureInfo["database"]]
                    } else if (match.signature.signatureLibraryRelease.library == "PIRSR") {
                        match.signature.signatureLibraryRelease.version = databaseInfo["PIRSR"]
                    }

                    // Handle PANTHER data
                    if (match.signature.signatureLibraryRelease.library == "PANTHER") {
                        updatePantherData(match, db_releases.panther.dirpath, panther_paint_dir, signatureAcc, entries)
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

def loadXRefFiles(prefix) {
    def iprFilePath = new File("${prefix}.ipr.json")
    def infoFilePath = new File("${prefix}.json")
    return [
        new ObjectMapper().readValue(iprFilePath, Map),
        new ObjectMapper().readValue(infoFilePath, Map)
    ]
}

def updatePantherData(def match, def pantherDir, def paintAnnoDir, def signatureAcc, def entries) {
    File paintAnnotationFile = new File("${pantherDir.toString()}/${paintAnnoDir}/${signatureAcc}.json")
    assert paintAnnotationFile.exists()
    def paintAnnotationsContent = new ObjectMapper().readValue(paintAnnotationFile, Map)
    String nodeId = match.treegrafter.ancestralNodeID
    def nodeData = paintAnnotationsContent[nodeId]
    if (nodeData != null) {
        match.treegrafter.subfamilyAccession = (nodeData[0] == "null") ? null : nodeData[0]
        match.treegrafter.proteinClass = (nodeData[2] == "null") ? null : nodeData[2]
        match.treegrafter.graftPoint = (nodeData[3] == "null") ? null : nodeData[3]
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
