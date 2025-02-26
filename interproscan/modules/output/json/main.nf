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
            jsonWriter.writeStartArray()  // start of results [...
            matchFiles.each { matchFile ->
                JsonReader.streamJson(matchFile.toString(), jacksonMapper) { String proteinMd5, JsonNode matchesNode ->
                    writeProtein(proteinMd5, matchesNode, proteinSeqData jsonWriter)
                }
            }
            jsonWriter.writeEndArray() // end of "results" ...]
        }
    }
}

def groupNucleotides(matchesFiles, seqDbPath) {
    /* Gather nucleotide Seq IDs and child protein MD5s so we can gather all ORFs from the same
    parent NT seq together in the final output. */
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

def writeProtein(String proteinMd5, JsonNode matchesNode, JsonGenerator jsonWriter) {
    /* Write data for a query protein sequence and its matches:
    { "sequence": sequence, "md5": proteinMd5, "matches": [], "xrefs": []}
    There may be multiple seqIds and desc for the same sequence/md5, use the first entry to get the seq. */

    // 1. {"sequence": seq, "md5": proteinMd5}
    proteinSeqData = getProteinSeqData(seqDbPath, proteinMd5)
    def sequence = proteinSeqData[0].split('\t')[2]
    jsonWriter.writeFieldName("sequence")
    jsonWriter.writeObject(sequence)
    jsonWriter.writeFieldName("md5")
    jsonWriter.writeObject(proteinMd5)

    // 2. {..., "matches": [{match}, {match}, {match}]}
    jsonWriter.writeFieldName("matches")
    jsonWriter.writeStartArray()
    matchesNode.fields().each { Map.Entry<String, JsonNode> entry ->
        JsonNode match = entry.value
        writeMatch(proteinMd5, match, jsonWriter)
    }
    jsonWriter.writeEndArray()

    // 3. {..., "xref": [{xref}, {xref}, {xref}]}
    writeProteinXref(proteinSeqData, jsonWriter)
}

def getProteinSeqData(def seqDbPath, String proteinMd5) {
    // Retrieve all associated seq IDs and desc for the given protein seq (id'd by it's md5 hash)
    try {
        def query = """
            SELECT P.id, P.description, S.sequence
            FROM PROTEIN AS P
            LEFT JOIN PROTEIN_SEQUENCE AS S ON P.protein_md5 = S.protein_md5
            WHERE P.protein_md5 = '$proteinMd5';"""
        def cmd = ["sqlite3", "--tabs", seqDbPath, query]
        def process = cmd.execute()
        process.waitFor()
        def output = process.in.text.trim()  // stndout
        return output.split("\\n")  // a protein seq may be associated with multiple seq IDs
    } catch (Exception e) {
        throw new Exception("Error when retrieving seq data for protein $proteinMd5: $e -- ${e.getCause()}")
    }
}

def writeMatch(String proteinMd5, JsonNode match, JsonGenerator jsonWriter) {
    // Write out an individual match to an array of matches. The structure is dependent on the memberDB.
    String memberDB = match.get("signature").get("signatureLibraryRelease").get("library").toLowerCase() ?: ""
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
            [
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
            [
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
            [
                "start"             : loc.get("start"),
                "end"               : loc.get("end"),
                "representative"    : loc.get("representative"),
                "evalue"            : loc.get("evalue"),
                "score"             : loc.get("score"),
                "location-fragments": loc.get("fragments"),
                "sites"             : loc.get("sites").collect { site ->
                    [
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
            [
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
            [
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
            [
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
            [
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
            [
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
            [
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
            [
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
            [
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
            [
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
            [
                "start"             : loc.get("start"),
                "end"               : loc.get("end"),
                "representative"    : loc.get("representative"),
                "hmmLength"         : loc.get("hmmLength"),
                "location-fragments": loc.get("fragments"),
            ]
        }
    ])
}

def writeProteinXref(def proteinSeqData, JsonGenerator jsonWriter) {
    /* "xref" : [ {
        "name" : "tr|A0A011PH51|A0A011PH51_9PROT OX=1454000",
        "id" : "tr|A0A011PH51|A0A011PH51_9PROT"
    } ] */
    jsonWriter.writeFieldName("xref")
    jsonWriter.writeStartArray()
    proteinSeqData.each { row ->
        row = row.split('\t')
        seqId = row[0]
        seqDesc = row[1]
        jsonWriter.writeObject(["name": "$seqId $seqDesc", "id": seqId])
    }
    jsonWriter.writeEndArray()
}


