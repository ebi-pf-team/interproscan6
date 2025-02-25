import com.fasterxml.jackson.core.JsonGenerator
import com.fasterxml.jackson.databind.JsonNode
import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.databind.SerializationFeature
import com.fasterxml.jackson.databind.node.ObjectNode
import java.util.regex.Pattern

process WRITE_JSON_OUTPUT {
    label 'local'

    input:
    val matchesFiles
    val outputPath
    val seqDbPath
    val nucleic
    val ips6Version

    exec:
    ObjectMapper jacksonMapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT)
    def NT_SEQ_ID_PATTERN = Pattern.compile(/^orf\d+\s+source=([^"]+)\s+coords=(\d+)\.\.(\d+)\s+length=(\d+)\s+frame=(\d+)\s+desc=(.*)$/)

    if (nucleic) {
        nucleicRelationships = groupNucleotides(matchesFiles, seqDbPath)
    } else {
        JsonWriter.streamJson(outputFilePath.toString(), jacksonMapper) { JsonGenerator jsonWriter ->
            // {"interproscan-version": str, "results": []}
            jsonWriter.writeStringField("interproscan-version", ips6Version)
            jsonWriter.writeFieldName("results")
            jsonWriter.writeStartArray()

            matchFiles.each { matchFile ->
                JsonReader.streamJson(matchFile.toString(), jacksonMapper) { String proteinMd5, JsonNode matchesNode ->
                    matchesNode.fields().each { Map.Entry<String, JsonNode> entry ->
                        JsonNode match = entry.value
                        writeMatch(proteinMd5, match, jsonWriter)
                    }
                }
            }
        }
    }

def writeMatch(String proteinMd5, JsonNode match, JsonGenerator jsonWriter) {
    /* Writes out a protein match = {
        "sequence": protein seq
        "md5": protein seq md5
        "matches": [ {match}, {match}, {match} ]
    } */
    String memberDB = match.get("signature").get("signatureLibraryRelease").get("library").toLowerCase() ?: ""
    // Write sequence
    // Write md5
    // add "matches" --> array
    switch (memberDB) {
        case "antifam":
            writeDefault(match, jsonWriter)
            break
        case "cath-gene3d":
            writeDefault(match, jsonWriter)
            break
        case "cath-funfam" || "funfam":
            writeDefault(match, jsonWriter)
            break
        case "cdd":
            writeCDD(match, jsonWriter)
            break
        case "coils":
            writeMinimalist(match, jsonWriter)
            break
        case "hamap":
            writeHAMAP(match, jsonWriter)
            break
        case "mobidb-lite" || "mobidb_lite":
            writeMobiDBlite(match, jsonWriter)
            break
        case "panther":
            writePANTHER(match, jsonWriter)
            break
        case "pfam":
            writeDefault(match, jsonWriter)
            break
        case "phobius":
            writeMinimalist(match, jsonWriter)
            break
        case "pirsf":
            writeDefault(match, jsonWriter)
            break
        case "pirsr":
            writeDefaultNoHmmBounds(match, jsonWriter)
            break
        case "prints":
            writePRINTS(match, jsonWriter)
            break
        case "prosite patterns":
            writePROSITEpatterns(match, jsonWriter)
            break
        case "prosite profiles":
            writePROSITEprofiles(match, jsonWriter)
            break
        case "signalp":
            writeSignalp(match, jsonWriter)
            break
        case "smart":
            writeSMART(match, jsonWriter)
            break
        case "superfamily":
            writeSUPERFAMILY(match, jsonWriter)
            break
        default:
            throw new UnsupportedOperationException("Unknown database '${memberDB}' for query protein with MD5 ${proteinMd5}")
    }
}

def writeDefault(JsonNode match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.get("signature"),
        "model-ac" : match.get("modelAccession"),
        "evalue"   : match.get("evalue"),
        "score"    : match.get("score"),
        "locations": match.get("locations").collect { loc ->
            return [
                "start"             : loc.get("start"),
                "end"               : loc.get("end"),
                "representative"    : loc.get("representative"),
                "hmmStart"          : loc.get("hmmStart"),
                "hmmEnd"            : loc.get("hmmEnd"),
                "hmmLength"         : loc.get("hmmLength"),
                "hmmBounds"         : loc.get("hmmBounds"),
                "evalue"            : loc.get("evalue"),
                "score"             : loc.get("score"),
                "envelopeStart"     : loc.get("envelopeStart"),
                "envelopeEnd"       : loc.get("envelopeEnd"),
                "location-fragments": loc.get("fragments"),
            ]
        }
    ])
}

def writeDefaultNoHmmBounds(JsonNode match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.get("signature"),
        "model-ac" : match.get("modelAccession"),
        "evalue"   : match.get("evalue"),
        "score"    : match.get("score"),
        "locations": match.get("locations").collect { loc ->
            return [
                "start"             : loc.get("start"),
                "end"               : loc.get("end"),
                "representative"    : loc.get("representative"),
                "hmmStart"          : loc.get("hmmStart"),
                "hmmEnd"            : loc.get("hmmEnd"),
                "hmmLength"         : loc.get("hmmLength"),
                "evalue"            : loc.get("evalue"),
                "score"             : loc.get("score"),
                "envelopeStart"     : loc.get("envelopeStart"),
                "envelopeEnd"       : loc.get("envelopeEnd"),
                "location-fragments": loc.get("fragments"),
            ]
        }
    ])
}

def writeCDD(JsonNode match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.get("signature"),
        "model-ac" : match.get("modelAccession"),
        "locations": match.get("locations").collect { loc ->
            return [
                "start"             : loc.get("start"),
                "end"               : loc.get("end"),
                "representative"    : loc.get("representative"),
                "evalue"            : loc.get("evalue"),
                "score"             : loc.get("score"),
                "location-fragments": loc.get("fragments"),
                "sites"             : loc.get("sites").collect { site ->
                    return [
                        "description"  : site.get("description"),
                        "numLocations" : site.get("numLocations"),
                        "siteLocations": site.get("siteLocations")
                    ]
                }
            ]
        }
    ])
}

def writeMinimalist(JsonNode match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.get("signature"),
        "model-ac" : match.get("modelAccession"),
        "locations": match.get("locations").collect { loc ->
            return [
                "start"             : loc.get("start"),
                "end"               : loc.get("end"),
                "representative"    : loc.get("representative"),
                "location-fragments": loc.get("fragments"),
            ]
        }
    ])
}

def writeHAMAP(JsonNode match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.get("signature"),
        "model-ac" : match.get("modelAccession"),
        "locations": match.get("locations").collect { loc ->
            return [
                "start"             : loc.get("start"),
                "end"               : loc.get("end"),
                "representative"    : loc.get("representative"),
                "score"             : loc.get("score"),
                "alignment"         : loc.get("targetAlignment"),
                "location-fragments": loc.get("fragments"),
            ]
        }
    ])
}

def writeMobiDBlite(JsonNode match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.get("signature"),
        "model-ac" : match.get("modelAccession"),
        "locations": match.get("locations").collect { loc ->
            return [
                "start"             : loc.get("start"),
                "end"               : loc.get("end"),
                "representative"    : loc.get("representative"),
                "location-fragments": loc.get("fragments"),
                "sequence-feature"  : loc.get("sequenceFeature"),
            ]
        }
    ])
}

def writePANTHER(JsonNode match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature"   : match.get("signature"),
        "model-ac"    : match.get("modelAccession"),
        "evalue"      : match.get("evalue"),
        "score"       : match.get("score"),
        "proteinClass": match.get("treegrafter").get("proteinClass"),
        "graftPoint"  : match.get("treegrafter").get("graftPoint"),
        "locations": match.get("locations").collect { loc ->
            return [
                "start"             : loc.get("start"),
                "end"               : loc.get("end"),
                "representative"    : loc.get("representative"),
                "hmmStart"          : loc.get("hmmStart"),
                "hmmEnd"            : loc.get("hmmEnd"),
                "hmmLength"         : loc.get("hmmLength"),
                "hmmBounds"         : loc.get("hmmBounds"),
                "evalue"            : loc.get("evalue"),
                "score"             : loc.get("score"),
                "envelopeStart"     : loc.get("envelopeStart"),
                "envelopeEnd"       : loc.get("envelopeEnd"),
                "location-fragments": loc.get("fragments"),
            ]
        }
    ])
}

def writePRINTS(JsonNode match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.get("signature"),
        "model-ac" : match.get("modelAccession"),
        "evalue"   : match.get("evalue"),
        "graphscan": match.get("graphscan"),
        "locations": match.get("locations").collect { loc ->
            return [
                "start"             : loc.get("start"),
                "end"               : loc.get("end"),
                "representative"    : loc.get("representative"),
                "pvalue"            : loc.get("pvalue"),
                "score"             : loc.get("score"),
                "motifNumber"       : loc.get("motifNumber"),
                "location-fragments": loc.get("fragments"),
            ]
        }
    ])
}

def writePROSITEpatterns(JsonNode match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.get("signature"),
        "model-ac" : match.get("modelAccession"),
        "locations": match.get("locations").collect { loc ->
            return [
                "start"             : loc.get("start"),
                "end"               : loc.get("end"),
                "representative"    : loc.get("representative"),
                "level"             : loc.get("level"),
                "cigarAlignment"    : loc.get("cigarAlignment"),
                "alignment"         : loc.get("targetAlignment"),
                "location-fragments": loc.get("fragments"),
            ]
        }
    ])
}

def writePROSITEprofiles(JsonNode match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.get("signature"),
        "model-ac" : match.get("modelAccession"),
        "locations": match.get("locations").collect { loc ->
            return [
                "start"             : loc.get("start"),
                "end"               : loc.get("end"),
                "representative"    : loc.get("representative"),
                "score"             : loc.get("score"),
                "alignment"         : loc.get("targetAlignment"),
                "location-fragments": loc.get("fragments"),
            ]
        }
    ])
}

def writeSignalp(JsonNode match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.get("signature"),
        "model-ac" : match.get("modelAccession"),
        "locations": match.get("locations").collect { loc ->
            return [
                "start"             : loc.get("start"),
                "end"               : loc.get("end"),
                "representative"    : loc.get("representative"),
                "pvalue"            : loc.get("pvalue"),
                "cleavageStart"     : loc.get("cleavageStart"),
                "cleavageEnd"       : loc.get("cleavageEnd"),
            ]
        }
    ])
}

def writeSMART(JsonNode match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.get("signature"),
        "model-ac" : match.get("modelAccession"),
        "evalue"   : match.get("evalue"),
        "score"     : match.get("score"),
        "locations": match.get("locations").collect { loc ->
            return [
                "start"             : loc.get("start"),
                "end"               : loc.get("end"),
                "representative"    : loc.get("representative"),
                "hmmStart"          : loc.get("hmmStart"),
                "hmmEnd"            : loc.get("hmmEnd"),
                "hmmLength"         : loc.get("hmmLength"),
                "hmmBounds"         : loc.get("hmmBounds"),
                "evalue"            : loc.get("evalue"),
                "score"             : loc.get("score"),
                "location-fragments": loc.get("fragments"),
            ]
        }
    ])
}

def writeSUPERFAMILY(JsonNode match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.get("signature"),
        "model-ac" : match.get("modelAccession"),
        "evalue"   : match.get("evalue"),
        "locations": match.get("locations").collect { loc ->
            return [
                "start"             : loc.get("start"),
                "end"               : loc.get("end"),
                "representative"    : loc.get("representative"),
                "hmmLength"         : loc.get("hmmLength"),
                "location-fragments": loc.get("fragments"),
            ]
        }
    ])
}






    def outputFilePath = "${outputPath}.ips6.json"
    JsonWriter.streamJson(outputFilePath.toString(), jacksonMapper) { JsonGenerator jsonWriter ->
        // {"interproscan-version": str, "results": []}
        jsonWriter.writeStringField("interproscan-version", ips6Version)
        jsonWriter.writeFieldName("results")
        jsonWriter.writeStartArray()

        JsonReader.streamArray(matches.toString(), jacksonMapper) { ObjectNode seqNode ->
            List seqMatches = []
            seqNode.get("matches").fields().each { matchNode ->
                Match matchObj = Match.fromJsonNode((JsonNode) matchNode.value)
                String memberDB = matchObj.signature.signatureLibraryRelease.library ?: ""
                matchResult = [
                    "signature": matchObj.signature,
                    "locations": []
                ]
                if (matchObj.locations) {
                    // Locations
                    matchObj.locations.each { location ->
                        locationResult = [
                            "start": location.start,
                            "end": location.end,
                            "representative": location.representative,
                            "location-fragments": location.fragments
                        ]
                        hmmBounds = location.getHmmBounds(location.hmmBounds)
                        def fields = membersLocationFields.get(memberDB, otherMembersLocationFields)
                        fields.each { field ->
                            switch (field) {
                                case "targetAlignment":
                                    locationResult["alignment"] = location.targetAlignment
                                    break
                                case "cigarAlignment":
                                    locationResult["cigarAlignment"] = location.cigarAlignment
                                    break
                                case "cleavageStart":
                                    locationResult["cleavageStart"] = matchObj.signalp.cleavageSiteStart
                                    break
                                case "cleavageEnd":
                                    locationResult["cleavageEnd"] = matchObj.signalp.cleavageSiteEnd
                                    break
                                case "envelopeStart":
                                    locationResult["envelopeStart"] = location.envelopeStart
                                    break
                                case "envelopeEnd":
                                    locationResult["envelopeEnd"] = location.envelopeEnd
                                    break
                                case "evalue":
                                    locationResult["evalue"] = location.evalue
                                    break
                                case "evalue-match":
                                    locationResult["evalue"] = matchObj.evalue
                                    break
                                case "hmmStart":
                                    locationResult["hmmStart"] = location.hmmStart
                                    break
                                case "hmmEnd":
                                    locationResult["hmmEnd"] = location.hmmEnd
                                    break
                                case "hmmLength":
                                    locationResult["hmmLength"] = location.hmmLength
                                    break
                                case "hmmBounds":
                                    locationResult["hmmBounds"] = hmmBounds
                                    break
                                case "level":
                                    locationResult["level"] = location.level
                                    break
                                case "motifNumber":
                                    locationResult["motifNumber"] = location.motifNumber
                                    break
                                case "pvalue":
                                    locationResult["pvalue"] = location.pvalue
                                    break
                                case "score":
                                    locationResult["score"] = location.score
                                    break
                                case "score-match":
                                    locationResult["score"] = matchObj.score
                                    break
                                case "sequence-feature":
                                    locationResult["sequence-feature"] = location.sequenceFeature
                                    break
                                default:
                                    println "Warning: Unknown fields '${field}'"
                                    break
                            }
                        }
                        // Sites
                        if (memberDB in ["PIRSR", "SFLD"]) {
                            locationResult["sites"] = location.sites ?: []
                        }
                        if (memberDB == "CDD") {
                            locationResult["sites"] = location.sites?.collect { site ->
                                site.remove("label")
                                site.remove("group")
                                site.remove("hmmStart")
                                site.remove("hmmEnd")
                                return site
                            } ?: []
                        }
                        matchResult["locations"].add(locationResult)
                    }
                }

                if (matchResult["locations"]) {
                    // Match level fields
                    if (memberDB in ["AntiFam", "FunFam", "CATH-Gene3D", "NCBIfam", "PANTHER", "Pfam", "PIRSF", "PIRSR", "SFLD", "SMART", "TMHMM"]) {
                        matchResult["evalue"] = matchObj.evalue
                        matchResult["score"] = matchObj.score
                    }
                    matchResult["model-ac"] = matchObj.modelAccession
                    switch (memberDB) {
                        case "SFLD":
                            matchResult["scope"] = null
                            break
                        case "PANTHER":
                            matchResult["name"] = matchObj.treegrafter.subfamilyDescription
                            matchResult["accession"] = matchObj.treegrafter.subfamilyAccession
                            matchResult["model-ac"] = matchObj.treegrafter.subfamilyAccession
                            matchResult["goXRefs"] = matchObj.signature?.entry?.goXRefs ?: []
                            matchResult["proteinClass"] = matchObj.treegrafter.proteinClass
                            matchResult["graftPoint"] = matchObj.treegrafter.graftPoint
                            break
                        case "PRINTS":
                            matchResult["evalue"] = matchObj.evalue
                            matchResult["graphscan"] = matchObj.graphScan
                            break
                        case ["SignalP_Euk", "SignalP-Prok"]:
                            matchResult["orgType"] = matchObj.signalp.orgType
                    }
                    seqMatches.add(matchResult)
                }
            }
            seqNode.set("matches", jacksonMapper.valueToTree(seqMatches))

            if (nucleic) {
                def nucleicSeqMd5 = seqNode.get("translatedFrom").get(0).get("md5")  // nucleic sequence md5 - same for all ORFs
                if (!nucleicResults.containsKey(nucleicSeqMd5)) {
                    nucleicResults[nucleicSeqMd5] = [
                        sequence          : seqNode.get("translatedFrom").get(0).get("sequence").asText(),
                        md5               : nucleicSeqMd5,
                        crossReferences   : [],
                        openReadingFrames : []
                    ]
                    seqNode.get("translatedFrom").forEach { ntSeq ->
                        nucleicResults[nucleicSeqMd5].crossReferences << [
                            name: "${ntSeq.get('id').asText()} ${ntSeq.get('description').asText()}",
                            id  : ntSeq.get('id').asText()
                        ]
                    }
                }
                def data = seqNode.get("xref").get(0).get("name").asText()
                def ntMatch = NT_SEQ_ID_PATTERN.matcher(seqNode.get("xref").get(0).get("name").asText())
                assert ntMatch.matches()
                nucleicResults[nucleicSeqMd5].openReadingFrames << [
                    start   : ntMatch.group(2) as int,
                    end     : ntMatch.group(3) as int,
                    strand  : (ntMatch.group(4) as int) < 4 ? "SENSE" : "ANTISENSE",
                    protein : [
                        sequence : seqNode.get("sequence").asText(),
                        md5      : seqNode.get("md5").asText(),
                        matches  : seqNode.get("matches"),
                        xref     : seqNode.get("xref")
                    ]
                ]
            }
            if (nucleic) {
                nucleicResults.values().each { result ->
                    jsonWriter.writeObject(result)
                }
            } else {
                jsonWriter.writeObject(seqNode)
            }
        }  // end of json reader
        jsonWriter.writeEndArray()
    }  // end of json writer
}

def groupNucleotides(matchesFiles, seqDbPath) {
    /* Gather nucleotide Seq IDs and child protein MD5s so we can gather all ORFs from the same parent NT seq together
    in the final output. */
    def nucleicRelationships = [:]  // [ntSeqId: [batchFilePath: [proteinMd5]]
    matchesFiles.each { matchFile ->
        JsonReader.streamJson(matchFile.toString(), jacksonMapper) { String proteinMd5, JsonNode matchesNode ->
            // get all parent NT seq Ids
            def query = """SELECT N.id
            FROM NUCLEOTIDE AS N
            LEFT JOIN PROTEIN_TO_NUCLEOTIDE AS N2P ON N.nt_md5 = N2P.nt_md5
            WHERE N2P.protein_md5 = '$proteinMd5'
            """
            def cmd = ["sqlite3", seqDbPath, query]
            def process = cmd.execute()
            process.waitFor()
            def output = process.in.text.trim() // stndout
            output = output.split("\\n")
            output.each { ntSeqId ->
                ntEntry = nucleicRelationships.computeIfAbsent(ntSeqId, { [:] } )
                relatedProteins = ntEntry.computeIfAbsent(matchFile.toString(), { [] as Set } )
                relatedProteins.add(proteinMd5)
            }
        }
    }
}
