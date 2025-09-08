import com.fasterxml.jackson.databind.ObjectMapper
import java.time.format.DateTimeFormatter
import java.time.LocalDate

process WRITE_TSV {
    label    'tiny'
    executor 'local'

    input:
    val matchesFiles
    val output_file
    val seqDbPath
    val nucleic

    exec:
    // Make connection to sequence Db through the JAVA interface for performance
    def db = new SeqDBQuery(seqDbPath.toString())

    // Prepare the output file
    def tsvFile = new File(output_file)
    tsvFile.text = "" // clear the file if it already exists

    // Set up buffer writing for faster output
    BATCH_SIZE = 10000
    lineBuffer = []
    def flushBuffer = {
        if (lineBuffer.size() > 0) {
            tsvFile.append(lineBuffer.join("\n") + "\n")
            lineBuffer.clear()
        }
    }

    def timingFile = new File("${output_file}_timing.log")
    timingFile.text = "timestamp,operation,duration_ms,details\n"

    // Each line contains: seqId md5 seqLength memberDb modelAcc sigDesc start end evalue status date entryAcc entryDesc goterms pathways
    def currentDate = LocalDate.now().format(DateTimeFormatter.ofPattern("dd-MM-yyyy"))
    Set<String> seenNucleicMd5s = new HashSet<>()

    matchesFiles.each { matchFile ->
        def jsonLoadStartTime = System.currentTimeMillis()
        Map proteins = new ObjectMapper().readValue(new File(matchFile.toString()), Map)
        def jsonLoadEndTime = System.currentTimeMillis()
        def jsonLoadDuration = jsonLoadEndTime - jsonLoadStartTime

        def timestamp = new Date().format("yyyy-MM-dd HH:mm:ss.SSS")
        def fileSize = new File(matchFile.toString()).length()
        timingFile.append("${timestamp},json_load,${jsonLoadDuration},\"${new File(matchFile.toString()).name} (${fileSize} bytes, ${proteins.size()} proteins)\"\n")

        if (nucleic) {
            nucleicToProteinMd5 = db.groupProteins(proteins)
            nucleicToProteinMd5.each { String nucleicMd5, Set<String> proteinMd5s ->
                if (!seenNucleicMd5s.contains(nucleicMd5)) {
                    seenNucleicMd5s.add(nucleicMd5)
                    def ntSeqData = db.nucleicMd5ToNucleicSeq(nucleicMd5)
                    ntSeqData.each { seq ->
                        String parentId = seq.id
                        proteinMd5s.each { String proteinMd5 ->
                            def proteinMatches = proteins[proteinMd5]
                            if (proteinMatches == null) return
                            proteinMatches.each { modelAcc, matchMap ->
                                def match = Match.fromMap(matchMap)
                                def proteinSeqData = db.getOrfSeq(proteinMd5, nucleicMd5)
                                writeMatch(proteinMd5, parentId, proteinSeqData, match, currentDate, tsvFile)
                            }
                        }
                    }   
                }
            }
        } else {
            def proteinMd5List = proteins.keySet().toList()
            
            // Time the database query
            def queryStartTime = System.currentTimeMillis()
            Map<String, List> seqData = db.proteinMd5sToProteinSeqs(proteinMd5List)
            def queryEndTime = System.currentTimeMillis()
            def queryDuration = queryEndTime - queryStartTime
            
            // Log query timing
            timestamp = new Date().format("yyyy-MM-dd HH:mm:ss.SSS")
            timingFile.append("${timestamp},db_query,${queryDuration},\"${proteinMd5List.size()} proteins from ${new File(matchFile.toString()).name}\"\n")
            
            // Time the entire proteins.each loop
            def proteinsLoopStartTime = System.currentTimeMillis()
            
            // Time protein processing
            // Don't change to `each`, for loops are faster and better optimised for the JVM
            for (Map.Entry entry : proteins.entrySet()) {
                String proteinMd5 = entry.key
                Map proteinMatches = entry.value
                def proteinStartTime = System.currentTimeMillis()
                
                proteinMatches.each { modelAcc, match ->
                    match = Match.fromMap(match)
                    writeMatch(proteinMd5, null, seqData[proteinMd5], match, currentDate,lineBuffer, BATCH_SIZE, flushBuffer)
                }
                
                def proteinEndTime = System.currentTimeMillis()
                def proteinDuration = proteinEndTime - proteinStartTime
                
                // Log protein processing timing
                timestamp = new Date().format("yyyy-MM-dd HH:mm:ss.SSS")
                timingFile.append("${timestamp},protein_processing,${proteinDuration},\"${proteinMd5} with ${proteinMatches.size()} matches\"\n")
            }
            
            def proteinsLoopEndTime = System.currentTimeMillis()
            def proteinsLoopDuration = proteinsLoopEndTime - proteinsLoopStartTime
            
            // Log total proteins.each loop timing
            timestamp = new Date().format("yyyy-MM-dd HH:mm:ss.SSS")
            timingFile.append("${timestamp},proteins_loop_total,${proteinsLoopDuration},\"${proteins.size()} proteins from ${new File(matchFile.toString()).name}\"\n")
        }
    }

    // Flush any remaining lines
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


// Keep this optimized Groovy formatLine method
def formatLine(String seqId, String seqMd5, int seqLength, Match match, Location loc, String currentDate) {
    // Extract data using native Groovy property access (fastest)
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
    
    // Optimized scoring value logic
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

    // Optimized GO terms building
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

    // Optimized pathways building
    def pathwayBuilder = new StringBuilder()
    boolean first = true
    for (pathway in interproPathways) {
        if (!first) pathwayBuilder.append('|')
        pathwayBuilder.append(pathway.databaseName).append(':').append(pathway.id)
        first = false
    }
    String pathwayString = pathwayBuilder.length() > 0 ? pathwayBuilder.toString() : "-"

    // Optimized final string building - pre-sized StringBuilder
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
