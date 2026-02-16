package uk.ac.ebi.interpro

import com.fasterxml.jackson.core.JsonFactory
import com.fasterxml.jackson.databind.ObjectMapper
import uk.ac.ebi.interpro.Entry
import uk.ac.ebi.interpro.GoXRefs
import uk.ac.ebi.interpro.Match
import uk.ac.ebi.interpro.PathwayXRefs
import uk.ac.ebi.interpro.RepresentativeInfo


class ProcessXrefs {
    static void run(String inputPath, Map dbReleases, boolean addGoTerms, 
                         boolean addPathways, String pantherPaintDirectory, String outputPath) {
        ObjectMapper mapper = new ObjectMapper()
        def (databaseInfo, entries, ipr2go, goInfo, ipr2pa, paInfo) = [null, null, null, null, null, null]
        if (dbReleases?.interpro?.dirpath != null) {
            String interproDir = dbReleases.interpro.dirpath.toString()
            File databasesJson = new File("${interproDir}/databases.json")
            databaseInfo = mapper.readValue(databasesJson, Map)
            File entriesJson = new File("${interproDir}/entries.json")
            if (entriesJson.exists()) {
                entries = mapper.readValue(entriesJson, Map)
            }
            if (addGoTerms) {
                (ipr2go, goInfo) = loadXRefFiles("${interproDir}/goterms")
            }
            if (addPathways) {
                (ipr2pa, paInfo) = loadXRefFiles("${interproDir}/pathways")
            }
        }

        def matches = mapper.readValue(new File(inputPath), Map.class)
        matches.each { String seqMd5, Map seqMatches ->
            seqMatches.each { modelAcc, matchMap ->
                Match match = Match.fromMap(matchMap)  // convert Map to Match object

                if (match.source != "InterPro-N") {
                    String signatureAcc = match.signature.accession
                    def signatureInfo = entries == null ? null : (entries[signatureAcc] ?: entries[modelAcc])

                    // Update library version
                    def version = (match.signature.signatureLibraryRelease.version == "null") ? null : match.signature.signatureLibraryRelease.version
                    if (!version && signatureInfo != null) {
                        match.signature.signatureLibraryRelease.version = databaseInfo?.get(signatureInfo["database"])
                    } else if (match.signature.signatureLibraryRelease.library == "PIRSR") {
                        match.signature.signatureLibraryRelease.version = databaseInfo?.get("PIRSR")
                    }

                    // Handle PANTHER data
                    if (match.signature.signatureLibraryRelease.library == "PANTHER") {
                        updatePantherData(match, dbReleases.panther.dirpath, pantherPaintDirectory, signatureAcc, entries, goInfo)
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

                    // update the match
                    seqMatches[modelAcc] = match
                }  // end of if (entries)
            } // end of seqMatches
        }  // end of Json reader / seq Id

        def jf = new JsonFactory()
        new File(outputPath).withWriter { writer ->
            def gen = jf.createGenerator(writer)
            mapper.writeValue(gen, matches)
            gen.close()
        }
    }

    static def loadXRefFiles(String prefix) {
        ObjectMapper mapper = new ObjectMapper()
        File iprFile = new File("${prefix}.ipr.json")
        File infoFile = new File("${prefix}.json")

        Map iprData = iprFile.exists() ? mapper.readValue(iprFile, Map) : null
        Map infoData = infoFile.exists() ? mapper.readValue(infoFile, Map) : null

        return [iprData, infoData]
    }

    static void updatePantherData(Match match, String pantherDir, String paintAnnoDir, 
                                 String signatureAcc, Map entries, Map goInfo) {
        Map<String,String> GO_PATTERN = ["P": "BIOLOGICAL_PROCESS", "C": "CELLULAR_COMPONENT", "F": "MOLECULAR_FUNCTION"]
        Map goTerms = goInfo == null ? null : goInfo["terms"]
        Map signatureEntry = entries == null ? null : entries[signatureAcc]
        File paintAnnotationFile = new File("${pantherDir.toString()}/${paintAnnoDir}/${signatureAcc}.json")
        assert paintAnnotationFile.exists()
        def paintAnnotationsContent = new ObjectMapper().readValue(paintAnnotationFile, Map)
        String nodeId = match.treegrafter.ancestralNodeID
        def nodeData = paintAnnotationsContent[nodeId]
        if (nodeData != null) {
            match.treegrafter.subfamilyAccession = (nodeData[0] == null) ? null : nodeData[0]
            if (nodeData[1] != null && goTerms != null) {
                nodeData[1].split(',').each { goTermId ->
                    def term = goTerms[goTermId]
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

        if (signatureEntry != null) {
            match.treegrafter.subfamilyName = signatureEntry["name"]
            match.treegrafter.subfamilyDescription = signatureEntry["description"]
        }
    }

    static void addXRefs(Match match, String interproAcc, def ipr2go, def goInfo, def ipr2pa, def paInfo) {
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
}
