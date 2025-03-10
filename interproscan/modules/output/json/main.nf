import com.fasterxml.jackson.core.JsonGenerator
import com.fasterxml.jackson.databind.JsonNode
import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.databind.SerializationFeature
import com.fasterxml.jackson.databind.node.ObjectNode
import java.util.regex.Pattern

process WRITE_JSON_OUTPUT {
    label 'local', 'ips6_container'

    input:
    val matchesFiles  // {query prot seq md5: {model acc: match}}
    val outputPath
    val seqDbPath
    val nucleic
    val ips6Version

    exec:
    ObjectMapper jacksonMapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT)
    def outputFilePath = "${outputPath}.ips6.json"

    SeqDatabase seqDatabase = new SeqDatabase(seqDbPath.toString())

    JsonWriter.streamJson(outputFilePath.toString(), jacksonMapper) { JsonGenerator jsonWriter ->
        // {"interproscan-version": str, "results": []}
        jsonWriter.writeStringField("interproscan-version", ips6Version)
        jsonWriter.writeFieldName("results")
        jsonWriter.writeStartArray()  // start of results [...
        matchesFiles.each { matchFile ->
            if (nucleic) {  // input was nucleic acid sequence
                proteinMatches = JsonReader.load(matchFile.toString(), jacksonMapper)
                nucleicToProteinMd5 = seqDatabase.groupProteins(proteinMatches)
                nucleicToProteinMd5.each { String nucleicMd5, Set<String> proteinMd5s ->
                    writeNucleic(nucleicMd5, proteinMd5s, proteinMatches, jsonWriter, seqDatabase)
                }
            } else {  // input was protein sequences
                JsonReader.streamJson(matchFile.toString(), jacksonMapper) { String proteinMd5, JsonNode matchesNode ->
                    writeProtein(proteinMd5, matchesNode, jsonWriter, seqDatabase)
                }
            }
        }
        jsonWriter.writeEndArray() // end of "results" ...]
    }
}

def writeNucleic(String nucleicMd5, Set<String> proteinMd5s, Map proteinMatches, JsonGenerator jsonWriter, SeqDatabase seqDatabase) {
    /* Write data for an input nucleic acid sequence, and then the matches for its associated ORFs
    {"sequence: nt seq, "md5": nt md5,
    "crossReferences": [{ntSeqData}, {ntSeqData}],
    "openReadingFrames": [{protein}, {protein}, {protein}]}
    There may be multiple nt seq Ids associated with the same nt seq, use the first entry to get the seq. */
    jsonWriter.writeStartObject()

    // 1. {"sequence": seq, "md5": ntMd5}
    ntSeqData = seqDatabase.getSeqData(nucleicMd5, true)
    String sequence = ntSeqData[0].split('\t')[-1]
    jsonWriter.writeStringField("sequence", sequence)
    jsonWriter.writeStringField("md5", nucleicMd5)
    
    // 2. {..., "crossReferences": [{ntSeqXref}, {ntSeqXref}]}
    jsonWriter.writeFieldName("crossReferences")
    writeXref(ntSeqData, jsonWriter)

    // 3. {..., "openReadingFrames": [{protein}, {protein}]}
    jsonWriter.writeFieldName("openReadingFrames")
    writeOpenReadingFrames(nucleicMd5, proteinMd5s, proteinMatches, jsonWriter, seqDatabase)

    jsonWriter.writeEndObject()
}

def writeOpenReadingFrames(String nucleicMd5, Set<String> proteinMd5s, Map proteinMatches, JsonGenerator jsonWriter, SeqDatabase seqDatabase){
    def SOURCE_NT_PATTERN = Pattern.compile(/^source=[^"]+\s+coords=(\d+)\.\.(\d+)\s+length=\d+\s+frame=(\d+)\s+desc=.*$/)

    jsonWriter.writeStartArray()
    proteinMd5s.each { String proteinMd5 ->
        // a proteinSeq/Md5 may be associated with multiple nt md5s/seq, only pull the data where the nt md5/seq is relevant
        proteinSeqData = seqDatabase.getSeqData(proteinMd5, false, nucleicMd5)
        proteinSeqData.each { row ->
            proteinDesc = row.split("\t")[1]
            def proteinSource = SOURCE_NT_PATTERN.matcher(proteinDesc)
            assert proteinSource.matches()
            jsonWriter.writeStartObject()
            jsonWriter.writeNumberField("start", proteinSource.group(1) as int)
            jsonWriter.writeNumberField("end", proteinSource.group(2) as int)
            jsonWriter.writeStringField("strand", (proteinSource.group(3) as int) < 4 ? "SENSE" : "ANTISENSE")
            jsonWriter.writeFieldName("protein")
            writeProtein(proteinMd5, proteinMatches[proteinMd5], jsonWriter, seqDatabase)
            jsonWriter.writeEndObject()
        }
    }
    jsonWriter.writeEndArray()
}

def writeProtein(String proteinMd5, JsonNode matchesNode, JsonGenerator jsonWriter, SeqDatabase seqDatabase) {
    /* Write data for a query protein sequence and its matches:
    { "sequence": sequence, "md5": proteinMd5, "matches": [], "xrefs": []}
    There may be multiple seqIds and desc for the same sequence/md5, use the first entry to get the seq. */
    jsonWriter.writeStartObject()

    // 1. {"sequence": seq, "md5": proteinMd5}
    proteinSeqData = seqDatabase.getSeqData(proteinMd5, false)
    String sequence = proteinSeqData[0].split('\t')[-1]
    jsonWriter.writeStringField("sequence", sequence)
    jsonWriter.writeStringField("md5", proteinMd5)

    // 2. {..., "matches": [{match}, {match}, {match}]}
    jsonWriter.writeFieldName("matches")
    jsonWriter.writeStartArray()
    matchesNode.fields().each { Map.Entry<String, JsonNode> entry ->
        JsonNode match = entry.value
        writeMatch(proteinMd5, match, jsonWriter)
    }
    jsonWriter.writeEndArray()

    // 3. {..., "xref": [{xref}, {xref}, {xref}]}
    jsonWriter.writeFieldName("xref")
    writeXref(proteinSeqData, jsonWriter)
    jsonWriter.writeEndObject()
}

def writeMatch(String proteinMd5, JsonNode match, JsonGenerator jsonWriter) {
    // Write out an individual match to an array of matches. The structure is dependent on the memberDB.
    String memberDB = JsonReader.asString(match.get("signature").get("signatureLibraryRelease").get("library")).toLowerCase() ?: ""
    switch (memberDB) {
        case "antifam":
            writeDefault(match, jsonWriter)
            break
        case "cath-gene3d":
            writeDefault(match, jsonWriter)
            break
        case "cath-funfam":
        case "funfam":  // use groovy case fall to allow multiple options
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
        case "mobidb-lite":
        case "mobidb_lite":  // use groovy case fall to allow multiple options
            writeMobiDBlite(match, jsonWriter)
            break
        case "ncbifam":
            writeDefault(match, jsonWriter)
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
        case "sfld":
            writeDefaultNoHmmBounds(match, jsonWriter)
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

def writeXref(seqData, JsonGenerator jsonWriter) {
    /* "xref"/"crossReferences" : [ {
        "name" : "tr|A0A011PH51|A0A011PH51_9PROT OX=1454000",
        "id" : "tr|A0A011PH51|A0A011PH51_9PROT"
    } ] */
    jsonWriter.writeStartArray()
    seqData.each { row ->
        row = row.split('\t')
        seqId = row[0]
        seqDesc = row[1]
        // jsonWrite.writeObject([name: "$seqId $seqDesc"]) does not correctly handle the formatted str
        jsonWriter.writeStartObject()
        jsonWriter.writeStringField("name", "$seqId $seqDesc")
        jsonWriter.writeStringField("id", seqId)
        jsonWriter.writeEndObject()
    }
    jsonWriter.writeEndArray()
}
