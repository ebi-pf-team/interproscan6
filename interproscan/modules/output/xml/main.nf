import com.fasterxml.jackson.databind.JsonNode
import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.databind.node.ObjectNode
import com.fasterxml.jackson.databind.SerializationFeature
import groovy.xml.MarkupBuilder
import java.time.format.DateTimeFormatter
import java.time.LocalDate
import java.io.StringWriter
import java.util.regex.Pattern

process WRITE_XML_OUTPUT {
    label 'local'

    input:
    val matches
    val outputPath
    val nucleic
    val ips6Version

    exec:
    def NT_SEQ_ID_PATTERN = Pattern.compile(/^orf\d+\s+?source=(.*?)?\s+coords=(\d+)\.\.(\d+)\s+length=\d+\s+frame=(\d+)\s+desc=(.*)$/)
    def writer = new StringWriter()
    def xml = new MarkupBuilder(writer)
    // set the correct encoding so symbols are formatted correctly in the final output
    xml.setEscapeAttributes(false) // Prevent escaping attributes
    xml.setEscapeText(false)       // Prevent escaping text

    ObjectMapper jacksonMapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT)
    // matches = [{sequence: str, md5: str, matches: {modelAcc: Match}, xref:[], translatedFrom: {}]
    String matchType = nucleic ? "nucleotide-sequence-matches" : "protein-matches"
    xml."$matchType"("interproscan-version": ips6Version) {
        if (nucleic) {
            /* Group by their parent NT seq first then write out to file.
            This will take a lot of memory but no longer will have the later standardisation of
            the sequence storage. Until then we can reduce memory requirements by storing data as Json/ObjectNodes */
            def processedNT = []
            def groupedByNT = [:]  // ntSeqMd5 : [ObjectNode protein]
            JsonReader.streamArray(matches.toString(), jacksonMapper) { ObjectNode seqNode ->
                // The seq node contains an translatedFrom for each identical nucleotide seq. The seqs are the same but id/desc differ
                if (seqNode.get("translatedFrom") != null) {
                    String ntSequenceMD5 = seqNode.get("translatedFrom").get(0).get("md5").asText()
                    if (groupedByNT.containsKey(ntSequenceMD5)) {
                        groupedByNT[ntSequenceMD5].add(seqNode)
                    } else {
                        groupedByNT[ntSequenceMD5] = [seqNode]
                    }
                }
            } // end of nucleotide JsonReader

            groupedByNT.each { String ntSequenceMD5, List<ObjectNode> seqNodes ->
                /* Multiple ORFs can be found in each input nucleotide sequence, and there can be multiple
                identical nucleotide sequences with different sequence ids. */
                if (!processedNT.contains(ntSequenceMD5)) {
                    processedNT << ntSequenceMD5
                    String ntSequence = seqNodes.get(0).get("translatedFrom").get(0).get("sequence").asText()
                    "nucleotide-sequence" {
                        sequence(md5: ntSequenceMD5, ntSequence)
                        // list all identical nucleotide seqs in the input file
                        seqNodes.get(0).get("translatedFrom").forEach { ntRef ->
                            xref(id: ntRef.get("id").asText(), name: "${ntRef.get('id').asText()} ${ntRef.get('description').asText().replaceAll('"', '')}")
                        }
                        seqNodes.forEach { ObjectNode proteinNode ->
                            def ntMatch = NT_SEQ_ID_PATTERN.matcher(proteinNode.get("xref").get(0).get("name").asText().replaceAll('"', ""))
                            assert ntMatch.matches()
                            def start = ntMatch.group(2) as int
                            def end = ntMatch.group(3) as int
                            def strand = (ntMatch.group(4) as int) < 4 ? "SENSE" : "ANTISENSE"
                            orf(start: start, end: end, strand: strand) {
                                protein {
                                    sequence(md5: proteinNode.get("md5").asText(), proteinNode.get("sequence").asText())
                                    matches {
                                        try {
                                            processMatches(proteinNode.get("matches"), xml)
                                        } catch (Exception e) {
                                            throw new Exception("Error processing XML:\n$e\n${e.printStackTrace()}\n${e.getCause()}", e)
                                        }
                                    }
                                    proteinNode.get("xref").forEach { ref ->
                                        xref(id: ref.get("id"), name: ref.get("name"))
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else {
            JsonReader.streamArray(matches.toString(), jacksonMapper) { ObjectNode seqNode ->
                protein {
                    try {
                        sequence(md5: seqNode.get("md5").asText(), seqNode.get("sequence").asText())
                        matches { processMatches(seqNode.get("matches"), xml) }
                        seqNode.get("xref").each { xrefData ->
                            xref(id: xrefData.get("id").asText(), name: xrefData.get("name").asText())
                        }
                    } catch (Exception e) {
                        println "Error processing XML:\n$e\n${e.printStackTrace()}\n${e.getCause()}"
                    }
                }
            }
        }
    }

    def outputFilePath = "${outputPath}.ips6.xml"
    new File(outputFilePath).text = writer.toString()
}

def processMatches(matches, xml) {
    List<String> hmmer3Members = ["AntiFam", "CATH-Gene3D", "FunFam", "hamap", "NCBIfam", "Pfam", "PIRSF", "PIRSR", "SFLD", "SUPERFAMILY"]
    List<String> hmmer3LocationFields = ["evalue", "score", "hmmStart", "hmmEnd", "hmmLength", "hmmBounds", "envelopeStart", "envelopeEnd"]
    Map<String, List<String>> memberLocationFields = [
        "CDD": ["match-evalue", "match-score"],
        "COILS": [],
        "HAMAP": ["score"],
        "MobiDB-lite": ["sequence-feature"],
        "PANTHER": ["hmmStart", "hmmEnd", "hmmLength", "hmmBounds", "envelopeStart", "envelopeEnd"],
        "Phobius": [],
        "PRINTS": ["motifNumber", "pvalue", "score"],
        "PROSITE patterns": ["level"],
        "PROSITE profiles": ["score"],
        "SignalP-Prok": ["score"],
        "SignalP-Euk": ["score"],
        "SMART": ["evalue", "score", "hmmStart", "hmmEnd", "hmmLength", "hmmBounds"],
        "SUPERFAMILY": ["evalue", "hmmLength"],
        "DeepTMHMM": []
    ]

    matches.fields().each { entry ->
        Match matchObj = Match.fromJsonNode(entry.value)
        def memberDb = matchObj.signature.signatureLibraryRelease.library
        def matchNodeName = hmmer3Members.contains(memberDb) ? "hmmer3" : (memberDb == "smart" ? "hmmer2" : memberDb)

        def matchAttributes = [:]
        if (hmmer3Members.findAll { it != "superfamily"}.contains(memberDb) || memberDb == "smart") {
            matchAttributes.evalue = matchObj.evalue
            matchAttributes.score = matchObj.score
        } else if (memberDb == "panther") {
            matchAttributes.ac = matchObj.treegrafter.subfamilyAccession
            matchAttributes.evalue = matchObj.evalue
            matchAttributes."graft-point" = matchObj.treegrafter.graftPoint
            matchAttributes.name = matchObj.signature.name
            matchAttributes.score = matchObj.score
        } else if (memberDb == "prints") {
            matchAttributes.evalue = matchObj.evalue
            matchAttributes.graftscan = matchObj.graphScan
        }

        xml."$matchNodeName-match"(matchAttributes) {
            def signatureAttributes = [ac: matchObj.signature.accession]
            def sigName = (matchObj.signature.name == "null") ? null : matchObj.signature.name
            def sigDesc = (matchObj.signature.description == "null") ? null : matchObj.signature.description
            if (sigName) { signatureAttributes.name = sigName }
            if (sigDesc) { signatureAttributes.desc = sigDesc }
            xml.signature(signatureAttributes) {
                xml.signatureLibraryRelease {
                    xml.library(matchObj.signature.signatureLibraryRelease.library)
                    xml.version(matchObj.signature.signatureLibraryRelease.version)
                }
                if (matchObj.signature.entry) {
                    Entry entryObj = matchObj.signature.entry
                    xml."entry"(
                        ac: entryObj.accession,
                        desc: (entryObj.description == null || entryObj.description == "null") ? "-" : entryObj.description,
                        name: (entryObj.name == null || entryObj.name == "null") ? "-" : entryObj.name,
                        type: (entryObj.type == null || entryObj.type == "null") ? "-" : entryObj.type,
                    ) {
                        if (!entryObj.goXRefs.isEmpty()) {
                            entryObj.goXRefs.each { goXrefObj ->
                                xml."go-xref"(
                                    category: goXrefObj.category,
                                    db: goXrefObj.databaseName,
                                    id: goXrefObj.id,
                                    name: goXrefObj.name
                                )
                            }
                        }
                        if (!entryObj.pathwayXRefs.isEmpty()) {
                            entryObj.pathwayXRefs.each { pathwayObj ->
                                xml."pathway-xref"(
                                    db: pathwayObj.databaseName,
                                    id: pathwayObj.id,
                                    name: pathwayObj.name
                                )
                            }
                        }
                    }
                }  // end of matchObj.signature.entry
            } // end of signature

            xml."model-ac"(memberDb == "panther" ? matchObj.treegrafter.subfamilyAccession : matchObj.modelAccession)

            if (matchObj.locations) {
                def fields = memberLocationFields.get(memberDb, hmmer3LocationFields)
                locations {
                    matchObj.locations.each { loc ->
                        def locationAttributes = getLocationAttributes(loc, fields, matchObj)
                        xml.location(locationAttributes) {
                            if (loc.fragments) {
                                "$matchNodeName-fragment" {
                                    loc.fragments.each { frag ->
                                        start(frag.start)
                                        end(frag.end)
                                        dcStatus(frag.dcStatus)
                                    }
                                }
                            }
                            if (memberDb in ["HAMAP", "PROSITE patterns", "PROSITE profiles"]) {
                                alignment(loc.targetAlignment ?: "")
                            }
                            if (loc.sites) {
                                "sites" {
                                    loc.sites.each { siteObj ->
                                        "$matchNodeName-site"(description: siteObj.description, numLocations: siteObj.numLocations) {
                                            if(siteObj.group){ group(siteObj.group) }
                                            if(siteObj.label){ label(siteObj.label) }
                                            if (memberDb != "cdd") {
                                                hmmStart(siteObj.hmmStart)
                                                hmmEnd(siteObj.hmmEnd)
                                            }
                                            xml."site-locations" {
                                                siteObj.siteLocations.each { siteLoc ->
                                                    xml."site-location"(
                                                        residue: siteLoc.residue,
                                                        start: siteLoc.start,
                                                        end: siteLoc.end
                                                    )
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

def getLocationAttributes(Location location, List<String> memberFields, Match matchObj) {
    def locationAttributes = [
        start: location.start,
        end: location.end,
        representative: location.representative
    ]
    memberFields.each { field ->
        switch (field) {
            case "evalue":
                locationAttributes.evalue = location.evalue
                break
            case "score":
                locationAttributes.score = location.score
                break
            case "match-evalue":
                locationAttributes.evalue = matchObj.evalue
                break
            case "match-score":
                locationAttributes.score = matchObj.score
                break
            case "motifNumber":
                locationAttributes.motifNumber = location.motifNumber
                break
            case "pvalue":
                locationAttributes.pvalue = location.pvalue
                break
            case "level":
                locationAttributes.level = location.level
                break
            case "sequence-feature":
                locationAttributes["sequence-feature"] = location.sequenceFeature
                break
            case "hmmStart":
                locationAttributes["hmm-start"] = location.hmmStart
                break
            case "hmmEnd":
                locationAttributes["hmm-end"] = location.hmmEnd
                break
            case "hmmLength":
                locationAttributes["hmm-length"] = location.hmmLength
                break
            case "hmmBounds":
                locationAttributes["hmm-bounds"] = location.getHmmBounds(location.hmmBounds)
                break
            case "envelopeStart":
                locationAttributes["env-start"] = location.envelopeStart
                break
            case "envelopeEnd":
                locationAttributes["env-end"] = location.envelopeEnd
                break
            default:
                println "Warning: Unknown fields '${field}'"
                break
        }
    }
    return locationAttributes
}