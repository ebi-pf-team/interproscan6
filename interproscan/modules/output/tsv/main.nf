import com.fasterxml.jackson.databind.JsonNode
import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.databind.SerializationFeature
import com.fasterxml.jackson.databind.node.ObjectNode
import java.time.format.DateTimeFormatter
import java.time.LocalDate

process WRITE_TSV_OUTPUT {
    label 'local'

    input:
    val matchesFiles
    val outputPath
    val nucleic

    exec:
    ObjectMapper jacksonMapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT)

    def tsvFile = new File("${outputPath}.ips6.tsv".toString())
    tsvFile.text = "" // clear the file if it already exists

    // Each line contains: seqId md5 seqLength memberDb modelAcc sigDesc start end evalue status date entryAcc entryDesc goterms pathways
    def currentDate = LocalDate.now().format(DateTimeFormatter.ofPattern("dd-MM-yyyy"))

    matchesFiles.each { matchFile ->
        JsonReader.streamJson(matchFile.toString(), jacksonMapper) { String proteinMd5, JsonNode matchNode ->

            Match match = Match.fromJsonNode((JsonNode) matchNode)
            String memberDb = match.signature.signatureLibraryRelease.library
            String sigDesc = (match.signature.description == null || match.signature.description == "null") ? '-' : match.signature.description
            String goterms = match.signature.entry?.goXRefs ?
                             (match.signature.entry.goXRefs.isEmpty() ? '-' :
                             match.signature.entry.goXRefs.collect { goXref -> "${goXref.id}(${goXref.databaseName})" }.join('|')) : '-'
            goterms = (goterms == "null") ? '-' : goterms
            String pathways = match.signature.entry?.pathwayXRefs ?
                             (match.signature.entry.pathwayXRefs.isEmpty() ? '-' :
                             match.signature.entry.pathwayXRefs.collect { ptXref -> "${ptXref.databaseName}:${ptXref.id}" }.join('|')) : '-'
            pathways = (pathways == "null") ? '-' : pathways
            String entryAcc = (match.signature.entry?.accession == null || match.signature.entry?.accession == "null") ? '-' : match.signature.entry?.accession
            String entryDesc = (match.signature.entry?.description == null || match.signature.entry?.description == "null") ? '-' : match.signature.entry?.description
            char status = 'T'

            seqIds = getSequenceData(proteinMd5, nucleic)
            seqIds.each { String seqId ->
                match.locations.each { Location loc ->
                    if ( nucleic ) {
                        seqNode.get("translatedFrom").forEach { translatedFromNode ->
                            seqId = "${translatedFromNode.get('id').asText()}_${xrefData.get('id').asText()}"
                            writeToTsv(tsvFile, match, currentDate, memberDb, sigDesc, goterms, pathways, entryAcc, entryDesc, status, seqNode, xrefData, loc, seqId)
                        }
                    } else {
                        writeToTsv(tsvFile, match, currentDate, memberDb, sigDesc, goterms, pathways, entryAcc, entryDesc, status, seqNode, xrefData, loc, seqId)
                    }
                }
            }
        }
    }
}

def getSequenceData(proteinMd5, nucleic) {
    if ( nucleic ) {
        return
    } else {
        return
    }
}

def writeToTsv(tsvFile, match, currentDate, memberDb,  sigDesc, goterms, pathways, entryAcc, entryDesc, status, seqNode, xrefData, loc, seqId) {
    int start = loc.start
    int end = loc.end
    def scoringValue = "-"
    switch (memberDb) {
        case ["CDD", "PRINT"]:
            scoringValue = match.evalue
            break
        case ["SignalP-Prok", "SignalP-Euk"]:
            scoringValue = loc.pvalue
            break
        case ["HAMAP", "PROSITE profiles"]:
            scoringValue = loc.score
            break
        case ["COILS", "MobiDB-lite", "Phobius", "PROSITE patterns", "DeepTMHMM"]:
            scoringValue = "-"
            break
        default:
            scoringValue = loc.evalue
            break
    }

    tsvFile.append([
        seqId, seqNode.get("md5").asText(), seqNode.get("sequence").asText().length(),
        memberDb,
        match.signature.accession, sigDesc,
        start, end, scoringValue, status,
        currentDate,
        entryAcc, entryDesc, goterms, "${pathways}\n"
    ].join('\t'))
}