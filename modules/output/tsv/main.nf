import com.fasterxml.jackson.databind.ObjectMapper
import java.time.format.DateTimeFormatter
import java.time.LocalDate

process WRITE_TSV {
    label    'mem_low', 'time_long'
    executor 'local'

    input:
    val matchesFiles
    val output_file
    val seqDbPath
    val nucleic

    exec:
    SeqDB db = new SeqDB(seqDbPath.toString())

    // Each line contains: seqId md5 seqLength memberDb modelAcc sigDesc start end evalue status date entryAcc entryDesc goterms pathways
    def currentDate = LocalDate.now().format(DateTimeFormatter.ofPattern("dd-MM-yyyy"))
    Set<String> seenNucleicMd5s = new HashSet<>()

    new File(output_file).withWriter { writer ->
        matchesFiles.each { matchFile ->
            Map proteins = new ObjectMapper().readValue(new File(matchFile.toString()), Map)
            if (nucleic) {
                def (nucleicToProteinMd5, ntSeqDataMap, orfDataMap) = 
                    db.retrieveAllNucleicSequenceData(proteins.keySet() as List)

                // Write TSV output
                nucleicToProteinMd5.each { String nucleicMd5, Set<String> proteinMd5s ->
                    if (!seenNucleicMd5s.contains(nucleicMd5)) {
                        seenNucleicMd5s.add(nucleicMd5)
                        def ntSeqData = ntSeqDataMap[nucleicMd5]
                        ntSeqData.each { seq ->
                            String parentId = seq.id
                            proteinMd5s.each { String proteinMd5 ->
                                def proteinMatches = proteins[proteinMd5]
                                if (proteinMatches == null) return

                                def proteinSeqData = orfDataMap[nucleicMd5][proteinMd5]
                                proteinMatches.each { modelAcc, matchMap ->
                                    def match = Match.fromMap(matchMap)
                                    writeMatch(proteinMd5, parentId, proteinSeqData, match, currentDate, writer)
                                }
                            }
                        }
                    }
                }
            } else {
                // Batch query all MD5s once
                def allMd5s = proteins.keySet() as List
                def md5ToSeqData = [:]

                allMd5s.collate(1000).each { batch ->
                    def result = db.proteinMd5ToProteinSeqs(batch)
                    md5ToSeqData.putAll(result)
                }

                // Now reuse the lookup map
                proteins.each { String proteinMd5, Map proteinMatches ->
                    def seqData = md5ToSeqData[proteinMd5]
                    proteinMatches.each { modelAcc, match ->
                        match = Match.fromMap(match)
                        writeMatch(proteinMd5, null, seqData, match, currentDate, writer)
                    }
                }
            }
        }
    }
}

def writeMatch(String proteinMd5, String proteinParentId, List seqData, Match match, String date, BufferedWriter writer) {
    seqData.each { row ->
        String seqId = proteinParentId ? "${proteinParentId}_${row.id}" : row.id
        int seqLength = row.sequence.trim().length()
        match.locations.each { Location loc ->
            String line = formatLine(seqId, proteinMd5, seqLength, match, loc, date)
            writer.writeLine(line)
        }
    }
}

def formatLine(String seqId, String seqMd5, int seqLength, Match match, Location loc, String currentDate) {
    String sigDesc = match.signature.description ?: '-'
    String entryAcc = match.signature.entry?.accession ?: '-'
    String entryDesc = match.signature.entry?.description ?: '-'
    def interproGoTerms = match.signature.entry?.goXRefs
    def interproPathways = match.signature.entry?.pathwayXRefs
    int start = loc.start
    int end = loc.end
    def scoringValue = "-"
    def pantherGoTerms = []
    String appl = (match.source == "InterPro-N" ? "InterPro-N" 
                                                : match.signature.signatureLibraryRelease.library
                                                ).toLowerCase()
    switch (appl) {
        case ["cdd", "prints"]:
            scoringValue = match.evalue
            break
        case ["signalp-prok", "signalp-euk"]:
            scoringValue = loc.pvalue
            break
        case ["hamap", "interpro-n", "prosite profiles"]:
            scoringValue = loc.score
            break
        case ["coils", "mobidb-lite", "phobius", "prosite patterns", "deeptmhmm", "tmbed"]:
            scoringValue = "-"
            break
        case "panther":
            pantherGoTerms = match.treegrafter.goXRefs.collect { "${it.id}(PANTHER)" }
            scoringValue = loc.evalue
            break
        default:
            scoringValue = loc.evalue
            break
    }

    goTerms = interproGoTerms.collect { "${it.id}(InterPro)" } + pantherGoTerms
    return [
        seqId, 
        seqMd5,
        seqLength,
        match.signature.signatureLibraryRelease.library,
        match.signature.accession,
        sigDesc,
        start,
        end,
        scoringValue,
        match.source,
        currentDate,
        entryAcc,
        entryDesc,
        goTerms.join("|") ?: "-",
        interproPathways.collect { "${it.databaseName}:${it.id}" }.join("|") ?: "-"
    ].join("\t")
}