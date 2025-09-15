import com.fasterxml.jackson.core.JsonGenerator
import com.fasterxml.jackson.databind.ObjectMapper
import groovy.json.JsonOutput

process XREFS {
    label    'tiny'
    executor 'local'

    input:
    tuple val(meta), val(matches_path)
    val db_releases
    val add_goterms
    val add_pathways
    val panther_paint_dir

    output:
    tuple val(meta), path("matches2xrefs.json")

    exec:
    def (databaseInfo, entries, ipr2go, goInfo, ipr2pa, paInfo) = [null, null, null, null, null, null]
    if (db_releases.interpro.dirpath != null) {
        String interproDir = db_releases.interpro.dirpath.toString()
        String databasesPath = "${interproDir}/databases.json"
        File databasesJson = new File(databasesPath)
        databaseInfo = new ObjectMapper().readValue(databasesJson, Map)
        String entriesPath = "${interproDir}/entries.json"
        File entriesJson = new File(entriesPath)
        entries = new ObjectMapper().readValue(entriesJson, Map)
        (ipr2go, goInfo) = loadXRefFiles("${interproDir}/goterms")
        if (add_pathways) {
            (ipr2pa, paInfo) = loadXRefFiles("${interproDir}/pathways")
        }
    }

    def matchesFileMap = new ObjectMapper().readValue(new File(matches_path.toString()), Map.class)
    matchesFileMap.each { String seqMd5, Map matches ->
        matches.each { modelAcc, matchMap ->
            match = Match.fromMap(matchMap)  // convert Map to Match object

            if (entries || entries != null) {
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
                    updatePantherData(match, db_releases.panther.dirpath, panther_paint_dir, signatureAcc, entries, add_goterms ? goInfo : null)
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
                         addXRefs(match, interproAcc, add_goterms ? ipr2go : null, goInfo, ipr2pa, paInfo)
                    }
                }

                // update the match
                matches[modelAcc] = match
            }  // end of if (entries)
        } // end of matches
    }  // end of Json reader / seq Id

    String outputFilePath = task.workDir.resolve("matches2xrefs.json")
    def json = JsonOutput.toJson(matchesFileMap)
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

def updatePantherData(def match, def pantherDir, def paintAnnoDir, def signatureAcc, def entries, def goInfo) {
    Map<String,String> GO_PATTERN = ["P": "BIOLOGICAL_PROCESS", "C": "CELLULAR_COMPONENT", "F": "MOLECULAR_FUNCTION"]
    File paintAnnotationFile = new File("${pantherDir.toString()}/${paintAnnoDir}/${signatureAcc}.json")
    assert paintAnnotationFile.exists()
    def paintAnnotationsContent = new ObjectMapper().readValue(paintAnnotationFile, Map)
    String nodeId = match.treegrafter.ancestralNodeID
    def nodeData = paintAnnotationsContent[nodeId]
    if (nodeData != null) {
        match.treegrafter.subfamilyAccession = (nodeData[0] == null) ? null : nodeData[0]
        if (nodeData[1] != null && goInfo?.terms != null) {
            nodeData[1].split(',').each { goTermId ->
                def term = goInfo.terms[goTermId]
                if (term) {
                    def (goTermName, goTermType) = term
                    match.treegrafter.addGoXRefs(
                        new GoXRefs(
                            goTermName,
                            "GO",
                            GO_PATTERN[goTermType],
                            goTermId
                        )
                    )
                }
            }
        }
        match.treegrafter.proteinClass = (nodeData[2] == null) ? null : nodeData[2]
        match.treegrafter.graftPoint = (nodeData[3] == null) ? null : nodeData[3]
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