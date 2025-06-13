import com.fasterxml.jackson.databind.ObjectMapper

import java.util.regex.Pattern

process WRITE_GFF3 {
    label    'tiny'
    executor 'local'

    input:
    val matchesFiles
    val outputPath
    val seqDbPath
    val nucleic
    val interproscan_version

    exec:
    SeqDB db = new SeqDB(seqDbPath.toString())
    def gff3File = new File("${outputPath}.gff3".toString())
    gff3File.text = "##gff-version 3.1.26\n"
    gff3File.append("##interproscan-version ${interproscan_version}\n")

    matchesFiles.each { matchFile ->
        matchFile = new File(matchFile.toString())
        Map proteins = new ObjectMapper().readValue(matchFile, Map)

        if (nucleic) {
            nucleicToProteinMd5 = db.groupProteins(proteins)
            nucleicToProteinMd5.each { String nucleicMd5, Set<String> proteinMd5s ->
                seqData = db.nucleicMd5ToNucleicSeq(nucleicMd5)
                seqId = seqData[0].id
                int seqLength = seqData[0].sequence.trim().length()
                gff3File.append("##sequence-region ${seqData[0].id} 1 ${seqLength}\n")

                proteinMd5s.each { String proteinMd5 ->
                    // a proteinSeq/Md5 may be associated with multiple nt md5s/seq, only pull the data where the nt md5/seq is relevant
                    proteinSeqData = db.getOrfSeq(proteinMd5, nucleicMd5)
                    proteinSeqData.each { row ->
                        gff3File.append(nucleotideLineFormat(seqId, row) + "\n")
                    }
                }
            }
        } else {
            proteins.each { String proteinMd5, Map matchesMap ->
                seqData = db.proteinMd5ToProteinSeq(proteinMd5)
                int seqLength = seqData[0].sequence.trim().length()
                gff3File.append("##sequence-region ${seqData[0].id} 1 ${seqLength}\n")

                matchesMap.each { modelAcc, match ->
                    match = Match.fromMap(match)
                    String memberDb = match.signature.signatureLibraryRelease.library
                    seqData.each { row ->
                        match.locations.each { Location loc ->
                            gff3File.append(proteinFormatLine(row, match, loc) + "\n")
                        }
                    }
                } // end of matches in matchesNode
            } // end of proteins.each
        } // end of nucleic else
    } // end of matchesFiles
}

def proteinFormatLine(seqInfo, match, loc) {

    String memberDb = match.signature.signatureLibraryRelease.library
    String entryAcc = match.signature.entry?.accession ?: '-'
    def goTerms = null
    if(memberDb == "PANTHER"){
        goTerms = match.treegrafter.goXRefs
    } else {
        goTerms = match.signature.entry?.goXRefs
    }

    def feature_type = null
    switch (memberDb) {
        case ["CATH-Gene3D", "CATH-FunFam", "CDD", "PROSITE profiles", "SMART", "SUPERFAMILY"]:
            feature_type = "polypeptide_domain"
            break
        case ["HAMAP", "MobiDB-lite", "PANTHER", "PIRSF", "PIRSR", "SFLD"]:
            feature_type = "polypeptide_region"
            break
        case ["NCBIFAM", "Pfam"]:
            feature_type = ["DOMAIN", "REPEAT"].contains(match.signature.type.toUpperCase()) ? "polypeptide_domain" : "polypeptide_region"
            break
        case "AntiFam":
            feature_type = "spurious_protein"
            break
        case "COILS":
            feature_type = "coiled_coil"
            break
        case "DeepTMHMM":
            feature_type = "transmembrane_helix"
            break
        case ["PRINTS", "PROSITE patterns"]:
            feature_type = "polypeptide_motif"
            break
        case ["SignalP-Prok", "SignalP-Euk"]:
            feature_type = "signal_peptide"
            break
        case "Phobius":
            feature_type = match.signature.type.toUpperCase() == "CYTOPLASMIC_DOMAIN" ? "cytoplasmic_polypeptide_region" :
                    match.signature.type.toUpperCase() == "NON_CYTOPLASMIC_DOMAIN" ? "non_cytoplasmic_polypeptide_region" :
                    match.signature.type.toUpperCase() == "TRANSMEMBRANE" ? "transmembrane_helix" :
                    "signal_peptide"
            break
    }

    def scoringValue = getScoringValue(match, loc)

    def attributes = [
            "Name=${match.signature.accession}",
            match.signature.description ? "Alias=${match.signature.description.replace(';',' ')}" : "Alias=${match.signature.name}",
            goTerms ? "Ontology_term=" + goTerms.collect { it.id }.join(",") : null,
            entryAcc && entryAcc != "-" ? "Dbxref=InterPro:${entryAcc}" : null,
            "type=${match.signature.type}",
            "representative=${loc.representative}",
    ].findAll { it }.join(";")

    return [
            seqInfo.id,
            memberDb,
            feature_type,
            loc.start,
            loc.end,
            scoringValue,
            ".", // strand
            ".", // phase
            attributes
    ].join("\t")
}

def nucleotideLineFormat(seqId, orfInfo){
    def SOURCE_NT_PATTERN = Pattern.compile(/^source=[^"]+\s+coords=(\d+)\.\.(\d+)\s+length=\d+\s+frame=(\d+)\s+desc=.*$/)
    def proteinSource = SOURCE_NT_PATTERN.matcher(orfInfo.description)
    assert proteinSource.matches()
    int start = proteinSource.group(1) as int
    int end = proteinSource.group(2) as int
    String strand = proteinSource.group(3) as int < 4 ? "+" : "-"

    return [
        seqId,
        "esl-translate",
        "CDS",
        strand == '-' ? end : start, // check if is reverse
        strand == '-' ? start : end, // check if is reverse
        ".", // score
        strand,
        0, //phase
        "ID=${seqId}_${orfInfo.id}"
    ].join("\t")
}

def getScoringValue(match, loc) {
    String memberDb = match.signature.signatureLibraryRelease.library

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
        case "PANTHER":
            scoringValue = loc.evalue
            break
        default:
            scoringValue = loc.evalue
    }
    def scoringValue = (scoringValue == "-" || scoringValue == null) ? "." : scoringValue
    return scoringValue
}