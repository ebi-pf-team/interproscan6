import com.fasterxml.jackson.databind.ObjectMapper

import java.util.regex.Pattern

process WRITE_GFF3 {
    label    'tiny'
    executor 'local'

    input:
    val matches_files
    val output_file
    val seq_db_file
    val nucleic
    val interproscan_version

    exec:
    def db = new SeqDBQuery(seq_db_file.toString())
    def gff3File = new File(output_file.toString())
    gff3File.text = "" // Clear file if it exists

    def tempFastaFile = new File("temp.fasta")
    tempFastaFile.text = ""

    def BATCH_SIZE = 5000
    def gff3LineBuffer = []
    def fastaLineBuffer = []

    def flushGff3Buffer = {
        if (gff3LineBuffer.size() > 0) {
            gff3File.append(gff3LineBuffer.join("\n") + "\n")
            gff3LineBuffer.clear()
        }
    }
    
    def flushFastaBuffer = {
        if (fastaLineBuffer.size() > 0) {
            tempFastaFile.append(fastaLineBuffer.join("\n") + "\n")
            fastaLineBuffer.clear()
        }
    }

    gff3File.append("##gff-version 3.1.26\n")
    gff3File.append("##interproscan-version ${interproscan_version}\n")

    matches_files.each { matchFile ->
        Map proteins = new ObjectMapper().readValue(new File(matchFile.toString()), Map)

        if (nucleic) {
            processNucleotides(db, proteins, gff3LineBuffer, fastaLineBuffer, BATCH_SIZE, flushGff3Buffer, flushFastaBuffer)
        } else {
            def proteinMd5List = proteins.keySet().toList()
            Map<String, List> seqData = db.proteinMd5sToProteinSeqs(proteinMd5List)

            for (Map.Entry entry : proteins.entrySet()) {
                String proteinMd5 = entry.key
                Map matchesMap = entry.value
                Map proteinSeqData = seqData[proteinMd5]

                String sequence = proteinSeqData[0].sequence.trim()
                int seqLength = sequence.length()

                proteinSeqData.each { row ->
                    gff3LineBuffer << "##sequence-region ${row.id} 1 ${seqLength}"

                    matchesMap.each { String modelAcc, Map match ->
                        Match match = Match.fromMap(matchMap)
                        for (Location loc : match.locations) {
                            String line = proteinFormatLine(row.id, match, loc, null, null, null)
                            gff3LineBuffer << line
                        }
                    }

                    if (gff3LineBuffer.size() >= BATCH_SIZE) {
                        flushGff3Buffer()
                        System.gc() // Make the jvm tidy up to prevent running out of memory
                    }

                    fastaLineBuffer << ">${row.id}"
                    fastaLineBuffer << "${sequence.replaceAll(/(.{60})/, '$1\n')}"

                    if (fastaLineBuffer.size() >= BATCH_SIZE) {
                        flushFastaBuffer()
                        System.gc() // Make the jvm tidy up to prevent running out of memory
                    }
                }         
            }
        }
    }

    flushGff3Buffer()
    flushFastaBuffer()

    gff3File.append("##FASTA\n")
    tempFastaFile.withReader { fastaReader ->
        fastaReader.eachLine { line ->
            gff3LineBuffer << "${line}\n"

            if (gff3LineBuffer.size() >= BATCH_SIZE) {
                flushGff3Buffer()
                System.gc() // Make the jvm tidy up to prevent running out of memory
            }
        }
    }

    flushGff3Buffer()

    tempFastaFile.delete()
    db.close()
}

def processNucleotidesBulkGFF3(SeqDBQuery db, Map proteins, Writer gff3Writer, Writer fastaWriter) {
    Pattern esl_pattern = Pattern.compile(/^source=[^"]+\s+coords=(\d+)\.\.(\d+)\s+length=\d+\s+frame=(\d+)\s+desc=.*$/)
    Set<String> seenNucleicMd5s = new HashSet<>()

    Set<String> allProteinMd5s = proteins.keySet().toSet()
    Map<String, Set<String>> nucleicToProteinMd5 = db.groupProteinsBulk(allProteinMd5s)
    
    // Filter to only unseen nucleotide MD5s
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

    newNucleicMd5s.each { String nucleicMd5 ->
        def ntSeqData = nucleotideSeqData[nucleicMd5]
        Set<String> proteinMd5sForNucleic = nucleicToProteinMd5[nucleicMd5]
        
        ntSeqData.each { seq ->
            String seqId = seq.id
            String sequence = seq.sequence.trim()
            int seqLength = sequence.length()
            gff3Writer.writeLine("##sequence-region ${seqId} 1 ${seqLength}")

            proteinMd5sForNucleic.each { String proteinMd5 ->
                def proteinMatches = proteins[proteinMd5]
                def proteinSeqData = orfSeqData[proteinMd5]?.get(nucleicMd5)
                
                proteinSeqData.each { row ->
                    def matcher = esl_pattern.matcher(row.description)
                    if (!matcher.matches()) return
                    
                    int start = matcher.group(1) as int
                    int end = matcher.group(2) as int
                    String strand = (matcher.group(3) as int) < 4 ? "+" : "-"
                    String parentId = "${seqId}_${row.id}"

                    String line
                    if (strand == "+") {
                        line = "${seqId}\tesl-translate\tCDS\t${start}\t${end}\t.\t${strand}\t0\tID=${parentId}"
                    } else {
                        line = "${seqId}\tesl-translate\tCDS\t${end}\t${start}\t.\t${strand}\t0\tID=${parentId}"
                    }

                    gff3Writer.writeLine(line)

                    proteinMatches.each { modelAcc, matchMap ->
                        def match = Match.fromMap(matchMap)
                        match.locations.each { Location loc ->
                            gff3Writer.writeLine(proteinFormatLine(seqId, match, loc, parentId, strand == "+" ? start : end, strand))
                        }
                    }
                }
            }

            fastaWriter.writeLine(">${seqId}")
            fastaWriter.writeLine("${sequence.replaceAll(/(.{60})/, '$1\n')}")  
        }
    }
}

def proteinFormatLine(seqId, match, loc, parentId, cdsStart, strand) {
    String memberDb = match.signature.signatureLibraryRelease.library

    def goTerms = []
    if(memberDb == "PANTHER" && match.treegrafter.goXRefs){
        goTerms += match.treegrafter.goXRefs
    }
    
    if (match.signature.entry?.goXRefs) {
        goTerms += match.signature.entry.goXRefs
    }

    def uniqueTermIds = new HashSet<String>()
    def goTermBuilder = new StringBuilder()
    boolean firstGo = true
    for (term in goTerms) {
        if (uniqueTermIds.add(term.id)) { // HashSet.add() returns false if already exists
            if (!firstGo) goTermBuilder.append(',')
            goTermBuilder.append(term.id)
            firstGo = false
        }
    }
    String goTermString = goTermBuilder.length() > 0 ? goTermBuilder.toString() : ""

    def feature_type = null
    switch (memberDb) {
        case ["CATH-Gene3D", "CATH-FunFam", "CDD", "PROSITE profiles", "SMART", "SUPERFAMILY"]:
            feature_type = "polypeptide_domain"
            break
        case ["NCBIFAM", "Pfam"]:
            feature_type = ["DOMAIN", "REPEAT"].contains(match.signature.type.toUpperCase()) ? "polypeptide_domain" : "polypeptide_region"
            break
        case ["PRINTS", "PROSITE patterns"]:
            feature_type = "polypeptide_motif"
            break
        case ["SignalP-Prok", "SignalP-Euk"]:
            feature_type = "signal_peptide"
            break
        case "AntiFam":
            feature_type = "spurious_protein"
            break
        case "COILS":
            feature_type = "coiled_coil"
            break
        case "TMHMM":
        case "DeepTMHMM":
            feature_type = match.signature.accession.toUpperCase() == "TRANSMEMBRANE ALPHA HELIX" ? "transmembrane_helix" :
                match.signature.accession.toUpperCase() == "TRANSMEMBRANE BETA BARREL" ? "transmembrane_polypeptide_region" :
                match.signature.accession.toUpperCase() == "PERIPLASMIC DOMAIN" ? "non_cytoplasmic_polypeptide_region" : "signal_peptide"
            break
        case "Phobius":
            feature_type = match.signature.name.toUpperCase() == "CYTOPLASMIC DOMAIN" ? "cytoplasmic_polypeptide_region" :
                match.signature.name.toUpperCase() == "NON CYTOPLASMIC DOMAIN" ? "non_cytoplasmic_polypeptide_region" :
                match.signature.name.toUpperCase() == "TRANSMEMBRANE REGION" ? "transmembrane_polypeptide_region" : "signal_peptide"
            break
        case "TMbed":
            feature_type = match.signature.accession.toUpperCase() == "TMHELIX_IN-TO-OUT" ? "transmembrane_helix" :
                match.signature.accession.toUpperCase() == "TMHELIX_OUT-TO-IN" ? "transmembrane_helix" :
                match.signature.accession.toUpperCase() == "TMBETA-OUT-TO-IN" ? "transmembrane_polypeptide_region" :
                match.signature.accession.toUpperCase() == "TMBETA-IN-TO-OUT" ? "transmembrane_polypeptide_region" : "signal_peptide"
            break
        default:
            // HAMAP, MobiDB-lite, Panther, PIRSF, PIRSR, SFLD
            feature_type = "polypeptide_region"
    }

    def score = null
    switch (memberDb) {
        case ["CDD", "PRINT"]:
            score = match.evalue
            break
        case ["HAMAP", "PROSITE profiles"]:
            score = loc.score
            break
        case ["SignalP-Prok", "SignalP-Euk"]:
            score = loc.pvalue
            break
        default:
            score = loc.evalue
    }

    String interproAccession = match.signature.entry?.accession

    def attributes = [
        match.signature.name ? "Name=${match.signature.name}" : null,
        "Alias=${match.signature.accession}",
        parentId ? "Parent=${parentId}" : null,
        interproAccession ? "Dbxref=InterPro:${interproAccession}" : null,
        goTermString ? "Ontology_term=${goTermString}" : null,
        "type=${match.signature.type}",
        "representative=${loc.representative}",
    ].findAll { it }

    StringBuilder sb = new StirngBuilder(256)
    sb.append(seqId).append("\t")
      .append(memberDb).append("\t")
      .append(feature_type).append("\t")
      .append(cdsStart ? (loc.start - 1) * 3 + cdsStart : loc.start).append("\t")
      .append(cdsStart ? loc.end * 3 + cdsStart -1 : loc.end).append("\t")
      .append(score ?: ".").append("\t")
      .append(strand ?: ".").append("\t")
      .append(parentId ? "0" : ".").append("\t")
      .append(attributes.join(";"))
    return sb.toString()
}
