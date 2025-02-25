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
    val seqDbPath
    val nucleic

    exec:
    ObjectMapper jacksonMapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT)

    def tsvFile = new File("${outputPath}.ips6.tsv".toString())
    tsvFile.text = "" // clear the file if it already exists

    // Each line contains: seqId md5 seqLength memberDb modelAcc sigDesc start end evalue status date entryAcc entryDesc goterms pathways
    def currentDate = LocalDate.now().format(DateTimeFormatter.ofPattern("dd-MM-yyyy"))

    matchesFiles.each { matchFile ->
        JsonReader.streamJson(matchFile.toString(), jacksonMapper) { String proteinMd5, JsonNode matchesNode ->
            matchesNode.fields().each { Map.Entry<String, JsonNode> entry ->
                Match match = Match.fromJsonNode((JsonNode) entry.value)
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

                seqData = getSeqData(seqDbPath, proteinMd5, nucleic)
                seqData.each { row ->  // Protein: [id, sequence]; Nucleic: [] from
                    match.locations.each { Location loc ->
                        if ( nucleic ) {
                            row = row.split('\t')
                            ntSeqId = row[0]
                            orfId = row[1]
                            seqLength = row[2].trim().length()
                            seqId = "${ntSeqId}_${orfId}"
                            writeToTsv(tsvFile, seqId, proteinMd5, seqLength, match, loc, memberDb, sigDesc, status, currentDate, entryAcc, entryDesc, goterms, pathways)
                        } else {
                            row = row.split('\t')
                            seqId = row[0]
                            seqLength = row[1].trim().length()
                            writeToTsv(tsvFile, seqId, proteinMd5, seqLength, match, loc, memberDb, sigDesc, status, currentDate, entryAcc, entryDesc, goterms, pathways)
                        }
                    }
                }
            } // end of matches in matchesNode
        } // end of JsonReader
    } // end of matchesFiles
}

def getSeqData(def seqDbPath, String querySeqMd5, boolean nucleic) {  // retrieve all seqIds associated with the query protein seq MD5 hash
    try {
        def query = ""
        if (nucleic) {
            query = """SELECT N.id, P.id, S.sequence
            FROM NUCLEOTIDE AS N
            LEFT JOIN PROTEIN_TO_NUCLEOTIDE AS N2P ON N.nt_md5 = N2P.nt_md5
            LEFT JOIN PROTEIN AS P ON N2P.protein_md5 = P.protein_md5
            LEFT JOIN PROTEIN_SEQUENCE AS S ON P.protein_md5 = N2P.protein_md5
            WHERE N2P.protein_md5 = '$querySeqMd5';"""
        } else {
            query = """SELECT P.id, S.sequence
            FROM PROTEIN AS P
            LEFT JOIN PROTEIN_SEQUENCE AS S ON P.protein_md5 = S.protein_md5
            WHERE P.protein_md5 = '$querySeqMd5';"""
        }
        def cmd = ["sqlite3", "--tabs", seqDbPath, query]
        def process = cmd.execute()
        process.waitFor()
        def output = process.in.text.trim()  // stndout
        return output.split("\\n")
    } catch (Exception e) {
        throw new Exception("Error when querying the internal seq DB: $e -- ${e.getCause()}", e)
    }
}

def writeToTsv(tsvFile, seqId, md5, seqLength, match, loc, memberDb, sigDesc, status, currentDate, entryAcc, entryDesc, goterms, pathways) {
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
        seqId, md5, seqLength, memberDb, match.signature.accession,
        sigDesc, start, end, scoringValue, status,
        currentDate, entryAcc, entryDesc, goterms, "${pathways}\n"
    ].join('\t'))
}