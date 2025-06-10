import com.fasterxml.jackson.databind.ObjectMapper

process WRITE_GFF3 {
    label    'tiny'
    executor 'local'

    input:
    val matchesFiles
    val outputPath
    val seqDbPath
    val nucleic

    exec:
    SeqDB db = new SeqDB(seqDbPath.toString())
    def gff3File = new File("${outputPath}.gff3".toString())
    gff3File.text = "##gff-version 3\n"

    matchesFiles.each { matchFile ->
        matchFile = new File(matchFile.toString())
        Map proteins = new ObjectMapper().readValue(matchFile, Map)
        proteins.each { String proteinMd5, Map matchesMap ->
            matchesMap.each { modelAcc, match ->
                match = Match.fromMap(match)
                String memberDb = match.signature.signatureLibraryRelease.library
                String sigDesc = match.signature.description ?: '-'
                def goterms = match.signature.entry?.goXRefs
                def pathways = match.signature.entry?.pathwayXRefs
                String entryAcc = match.signature.entry?.accession ?: '-'
                String entryName = match.signature.entry?.name ?: '-'
                String entryDesc = match.signature.entry?.description ?: '-'
                String entryType = match.signature.entry?.type ?: '-'
                seqData = nucleic ? db.proteinMd5ToNucleicSeq(proteinMd5) : db.proteinMd5ToProteinSeq(proteinMd5)
                seqData.each { row ->  // Protein or Nucleic: [id, desc, sequence]
                    String seqId = nucleic ? "${row.nid}_${row.pid}" : row.id
                    int seqLength = row.sequence.trim().length()
                    match.locations.each { Location loc ->
                        gff3File.append(formatLine(
                                seqId, match, loc, memberDb, sigDesc,
                                entryAcc, entryName, entryDesc, entryType,
                                goterms, pathways
                        ) + "\n")
                    }
                }
            } // end of matches in matchesNode
        } // end of proteins.each
    } // end of matchesFiles
}

def formatLine(
        seqId,
        match,
        loc,
        memberDb,
        sigDesc,
        entryAcc,
        entryName,
        entryDesc,
        entryType,
        interproGoTerms,
        interproPathways) {

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

    def scoringValue
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
    scoringValue = (scoringValue == "-" || scoringValue == null) ? "." : scoringValue

    def attributes = [
            "Name=${match.signature.accession}",
            entryDesc ? "Alias=${entryDesc}" : "Alias=${entryName}",
            interproGoTerms.collect { "Ontology_term=${it.id}" }.join(","),
            entryAcc ? "Dbxref=InterPro:${entryAcc}" : null,
            "type=${match.signature.type}",
            "representative=${loc.representative}",
    ].findAll { it }.join(";")

    return [
            seqId,
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
