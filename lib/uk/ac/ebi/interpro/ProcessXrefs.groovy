package uk.ac.ebi.interpro

import java.nio.file.Path
import com.fasterxml.jackson.core.JsonFactory
import com.fasterxml.jackson.databind.ObjectMapper
import uk.ac.ebi.interpro.Entry
import uk.ac.ebi.interpro.GoXRefs
import uk.ac.ebi.interpro.Match
import uk.ac.ebi.interpro.PathwayXRefs
import uk.ac.ebi.interpro.RepresentativeInfo


class ProcessXrefs {
    static void run(Path inputPath, Map dbReleases, boolean addGoTerms, 
                    boolean addPathways, String pantherPaintDirectory, Path outputPath) {
        ObjectMapper mapper = new ObjectMapper()
        def (databaseInfo, entries, ipr2go, goInfo, ipr2pa, paInfo) = [null, null, null, null, null, null]
        if (dbReleases?.interpro?.dirpath != null) {
            Path interproDir = dbReleases.interpro.dirpath
            Path databasesJson = interproDir.resolve("databases.json")
            databaseInfo = databasesJson.newReader().withCloseable { reader -> 
                mapper.readValue(reader, Map)
            }
            Path entriesJson = interproDir.resolve("entries.json")
            if (entriesJson.exists()) {
                entries = entriesJson.newReader().withCloseable { reader -> 
                    mapper.readValue(reader, Map)
                }
            }
            if (addGoTerms) {
                (ipr2go, goInfo) = loadXRefFiles(interproDir.resolve("goterms"))
            }
            if (addPathways) {
                (ipr2go, goInfo) = loadXRefFiles(interproDir.resolve("pathways"))
            }
        }

        def matches = inputPath.newReader().withCloseable { reader ->
            mapper.readValue(reader, Map.class)
        }
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
                        Path pantherDir = dbReleases.panther.dirpath
                        updatePantherData(match, pantherDir, pantherPaintDirectory, signatureAcc, entries, goInfo)
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
        outputPath.withWriter { writer ->
            def gen = jf.createGenerator(writer)
            mapper.writeValue(gen, matches)
            gen.close()
        }
    }

    static def loadXRefFiles(Path prefix) {
        Path iprFile = prefix.resolveSibling(prefix.fileName.toString() + ".ipr.json")
        Path infoFile = prefix.resolveSibling(prefix.fileName.toString() + ".json")
        ObjectMapper mapper = new ObjectMapper()
    
        Map iprData = null
        if (iprFile.exists()) {
            iprData = iprFile.newReader().withCloseable { reader ->
                mapper.readValue(reader, Map)
            }
        }

        Map infoData = null
        if (infoFile.exists()) {
            iprData = infoFile.newReader().withCloseable { reader ->
                mapper.readValue(reader, Map)
            }
        }

        return [iprData, infoData]
    }

    static void updatePantherData(Match match, Path pantherDir, String paintAnnoDir, 
                                 String signatureAcc, Map entries, Map goInfo) {
        Map<String,String> GO_PATTERN = ["P": "BIOLOGICAL_PROCESS", "C": "CELLULAR_COMPONENT", "F": "MOLECULAR_FUNCTION"]
        Map goTerms = goInfo == null ? null : goInfo["terms"]
        Map signatureEntry = entries == null ? null : entries[signatureAcc]

        Path paintAnnotationFile = pantherDir.resolve("${paintAnnoDir}/${signatureAcc}.json")
        assert paintAnnotationFile.exists()

        def paintAnnotationsContent = paintAnnotationFile.newReader().withCloseable { reader ->
            new ObjectMapper().readValue(reader, Map)
        }
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
