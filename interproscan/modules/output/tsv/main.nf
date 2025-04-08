import com.fasterxml.jackson.databind.ObjectMapper
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
    SeqDB db = new SeqDB(seqDbPath.toString())
    def tsvFile = new File("${outputPath}.ips6.tsv".toString())
    tsvFile.text = "" // clear the file if it already exists

    // Each line contains: seqId md5 seqLength memberDb modelAcc sigDesc start end evalue status date entryAcc entryDesc goterms pathways
    def currentDate = LocalDate.now().format(DateTimeFormatter.ofPattern("dd-MM-yyyy"))

    matchesFiles.each { matchFile ->
        matchFile = new File(matchFile.toString())
        Map proteins = new ObjectMapper().readValue(matchFile, Map)
        proteins.each { String proteinMd5, Map matchesMap ->
            matchesMap.each { modelAcc, match ->
                match = Match.fromMap(match)
                String memberDb = match.signature.signatureLibraryRelease.library
                String sigDesc = match.signature.description ?: '-'
                String goterms = match.signature.entry?.goXRefs ? match.signature.entry.goXRefs.collect { "${it.id}(${it.databaseName})" }.join('|') : '-'
                String pathways = match.signature.entry?.pathwayXRefs ? match.signature.entry.pathwayXRefs.collect { "${it.databaseName}:${it.id}" }.join('|') : '-'
                String entryAcc = match.signature.entry?.accession ?: '-'
                String entryDesc = match.signature.entry?.description ?: '-'
                char status = 'T'
                seqData = nucleic ? db.proteinMd5ToNucleicSeq(proteinMd5) : db.proteinMd5ToProteinSeq(proteinMd5)
                seqData.each { row ->  // Protein or Nucleic: [id, desc, sequence]
                    String seqId = nucleic ? "${row.nid}_${row.pid}" : row.id
                    int seqLength = row.sequence.trim().length()
                    match.locations.each { Location loc ->
                        writeToTsv(tsvFile, seqId, proteinMd5, seqLength, match, loc, memberDb, sigDesc, status, currentDate, entryAcc, entryDesc, goterms, pathways)
                    }
                }
            } // end of matches in matchesNode
        } // end of proteins.each
    } // end of matchesFiles
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