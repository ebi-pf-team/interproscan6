import com.fasterxml.jackson.databind.ObjectMapper
import java.time.format.DateTimeFormatter
import java.time.LocalDate

process WRITE_TSV {
    label    'tiny'
    executor 'local'

    input:
    val matchesFiles
    val output_file
    val seq_db_file
    val nucleic

    exec:
    def db = new SeqDBQuery(seq_db_file.toString())

    def tsvFile = new File(output_file)
    tsvFile.text = "" // clear the file if it already exists

    BATCH_SIZE = 5000
    lineBuffer = []
    def flushBuffer = {
        if (lineBuffer.size() > 0) {
            tsvFile.append(lineBuffer.join("\n") + "\n")
            lineBuffer.clear()
        }
    }

    // Each line contains: seqId md5 seqLength memberDb modelAcc sigDesc start end evalue status date entryAcc entryDesc goterms pathways
    def currentDate = LocalDate.now().format(DateTimeFormatter.ofPattern("dd-MM-yyyy"))
    Set<String> seenNucleicMd5s = new HashSet<>()

    matchesFiles.each { matchFile ->
        Map proteins = new ObjectMapper().readValue(new File(matchFile.toString()), Map)

        if (nucleic) {
            processNucleotidesBulk(javaDb, proteins, seenNucleicMd5s, currentDate, lineBuffer, BATCH_SIZE, flushBuffer)
        } else {
            def proteinMd5List = proteins.keySet().toList()
            Map<String, List> seqData = db.proteinMd5sToProteinSeqs(proteinMd5List)
    
            // Don't change to `each`, for loops are faster and better optimised for the JVM
            for (Map.Entry entry : proteins.entrySet()) {
                String proteinMd5 = entry.key
                Map proteinMatches = entry.value
                proteinMatches.each { modelAcc, match ->
                    match = Match.fromMap(match)
                    writeMatch(proteinMd5, null, seqData[proteinMd5], match, currentDate,lineBuffer, BATCH_SIZE, flushBuffer)
                }
            }
        }
    }

    flushBuffer()
}

def writeMatch(String proteinMd5, String proteinParentId, List seqData, Match match, String date, List lineBuffer, int batchSize, Closure flushBuffer) {
    seqData.each { seq ->
        String seqId = proteinParentId ? "${proteinParentId}_${seq.id}" : seq.id
        int seqLength = seq.sequence.trim().length()
        for (loc in match.locations) {
            String line = formatLine(seqId, proteinMd5, seqLength, match, loc, date)
            lineBuffer << line
            if (lineBuffer.size() >= batchSize) {
                flushBuffer()
            }
        }
    }
}


def formatLine(String seqId, String seqMd5, int seqLength, Match match, Location loc, String currentDate) {
    String memberDb = match.signature.signatureLibraryRelease.library
    String sigDesc = match.signature.description ?: '-'
    String accession = match.signature.accession
    String entryAcc = match.signature.entry?.accession ?: '-'
    String entryDesc = match.signature.entry?.description ?: '-'
    
    def interproGoTerms = match.signature.entry?.goXRefs ?: []
    def interproPathways = match.signature.entry?.pathwayXRefs ?: []
    int start = loc.start
    int end = loc.end
    def scoringValue = "-"
    def pantherGoTerms = []
    
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
        case ["COILS", "MobiDB-lite", "Phobius", "PROSITE patterns", "DeepTMHMM", "TMbed"]:
            scoringValue = "-"
            break
        case "PANTHER":
            if (match.treegrafter?.goXRefs) {
                pantherGoTerms = match.treegrafter.goXRefs.collect { "${it.id}(PANTHER)" }
            }
            scoringValue = loc.evalue
            break
        default:
            scoringValue = loc.evalue
            break
    }

    def goTermBuilder = new StringBuilder()
    boolean firstGo = true
    for (term in interproGoTerms) {
        if (!firstGo) goTermBuilder.append('|')
        goTermBuilder.append(term.id).append("(InterPro)")
        firstGo = false
    }
    for (term in pantherGoTerms) {
        if (!firstGo) goTermBuilder.append('|')
        goTermBuilder.append(term)
        firstGo = false
    }
    String goTermString = goTermBuilder.length() > 0 ? goTermBuilder.toString() : "-"

    def pathwayBuilder = new StringBuilder()
    boolean first = true
    for (pathway in interproPathways) {
        if (!first) pathwayBuilder.append('|')
        pathwayBuilder.append(pathway.databaseName).append(':').append(pathway.id)
        first = false
    }
    String pathwayString = pathwayBuilder.length() > 0 ? pathwayBuilder.toString() : "-"

    StringBuilder sb = new StringBuilder(256)
    sb.append(seqId).append('\t')
      .append(seqMd5).append('\t')
      .append(seqLength).append('\t')
      .append(memberDb).append('\t')
      .append(accession).append('\t')
      .append(sigDesc).append('\t')
      .append(start).append('\t')
      .append(end).append('\t')
      .append(scoringValue).append('\t')
      .append("T").append('\t')
      .append(currentDate).append('\t')
      .append(entryAcc).append('\t')
      .append(entryDesc).append('\t')
      .append(goTermString).append('\t')
      .append(pathwayString)
    return sb.toString()
}

def processNucleotidesBulk(SeqDBQuery db, Map proteins, Set seenNucleicMd5s, String currentDate, 
                          List lineBuffer, int batchSize, Closure flushBuffer) {
    
    Set<String> allProteinMd5s = proteins.keySet().toSet()
    Map<String, Set<String>> nucleicToProteinMd5 = db.groupProteinsBulk(allProteinMd5s)
    List<String> newNucleicMd5s = nucleicToProteinMd5.keySet().findAll { !seenNucleicMd5s.contains(it) }
    seenNucleicMd5s.addAll(newNucleicMd5s)
    if (newNucleicMd5s.isEmpty()) {
        return
    }

    Map<String, List> nucleotideSeqData = db.nucleicMd5sToNucleicSeqs(newNucleicMd5s)

    // Get all ORF sequences for relevant protein/nucleotide combinations
    List<String> relevantProteinMd5s = []
    newNucleicMd5s.each { nucleicMd5 ->
        relevantProteinMd5s.addAll(nucleicToProteinMd5[nucleicMd5])
    }
    Map<String, Map<String, List>> orfSeqData = db.getOrfSeqsBulk(relevantProteinMd5s, newNucleicMd5s)

    // Process the results
    newNucleicMd5s.each { String nucleicMd5 ->
        def ntSeqData = nucleotideSeqData[nucleicMd5]
        if (!ntSeqData) return
        
        Set<String> proteinMd5sForNucleic = nucleicToProteinMd5[nucleicMd5]
        
        ntSeqData.each { seq ->
            String parentId = seq.id
            proteinMd5sForNucleic.each { String proteinMd5 ->
                def proteinMatches = proteins[proteinMd5]
                if (proteinMatches == null) return
                
                proteinMatches.each { modelAcc, matchMap ->
                    def match = Match.fromMap(matchMap)
                    def proteinSeqData = orfSeqData[proteinMd5]?.get(nucleicMd5)
                    if (proteinSeqData) {
                        writeMatch(proteinMd5, parentId, proteinSeqData, match, currentDate, lineBuffer, batchSize, flushBuffer)
                    }
                }
            }
        }
    }
}
