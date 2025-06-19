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
    SeqDB db = new SeqDB(seq_db_file.toString())
    def gff3File = new File(output_file.toString())
    gff3File.text = "##gff-version 3.1.26\n"
    gff3File.append("##interproscan-version ${interproscan_version}\n")

    Pattern esl_pattern = Pattern.compile(/^source=[^"]+\s+coords=(\d+)\.\.(\d+)\s+length=\d+\s+frame=(\d+)\s+desc=.*$/)

    matches_files.each { matchFile ->
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
                        def matcher = esl_pattern.matcher(row.description)
                        assert matcher.matches()
                        int start = matcher.group(1) as int
                        int end = matcher.group(2) as int
                        String strand = (matcher.group(3) as int) < 4 ? "+" : "-"
                        String parentId = "${seqId}_${row.id}"

                        String line
                        if (strand == "+") {
                            line = "${seqId}\tesl-translate\tCDS\t${start}\t${end}\t.\t${strand}\t0\tID=${parentId}\n"
                        } else {
                            line = "${seqId}\tesl-translate\tCDS\t${end}\t${start}\t.\t${strand}\t0\tID=${parentId}\n"
                        }

                        gff3File.append(line)

                        proteins[proteinMd5].each { modelAcc, match->
                            match = Match.fromMap(match)
                            match.locations.each { Location loc ->
                                gff3File.append(proteinFormatLine(seqId, match, loc, parentId, strand == "+" ? start : end, strand) + "\n")
                            }
                        }
                    }
                }
            }
        } else {
            proteins.each { String proteinMd5, Map matchesMap ->
                seqData = db.proteinMd5ToProteinSeq(proteinMd5)
                int seqLength = seqData[0].sequence.trim().length()

                seqData.each { row ->
                    gff3File.append("##sequence-region ${row.id} 1 ${seqLength}\n")

                    matchesMap.each { modelAcc, match ->
                        match = Match.fromMap(match)
                        
                        match.locations.each { Location loc ->
                            gff3File.append(proteinFormatLine(row.id, match, loc, null, null, null) + "\n")
                        }
                    }
                }                
            }
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

    def uniqueTerms = [:]
    goTerms.each { term -> uniqueTerms[term.id] = term }
    goTerms = uniqueTerms.values() as List

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
        case "DeepTMHMM":
            feature_type = "transmembrane_helix"
            break
        case "Phobius":
            feature_type = match.signature.type.toUpperCase() == "CYTOPLASMIC_DOMAIN" ? "cytoplasmic_polypeptide_region" :
                    match.signature.type.toUpperCase() == "NON_CYTOPLASMIC_DOMAIN" ? "non_cytoplasmic_polypeptide_region" :
                    match.signature.type.toUpperCase() == "TRANSMEMBRANE" ? "transmembrane_helix" : "signal_peptide"
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
        goTerms ? "Ontology_term=${goTerms.collect{ it.id }.join(',')}" : null,
        "type=${match.signature.type}",
        "representative=${loc.representative}",
    ].findAll { it }

    return [
        seqId,
        memberDb,
        feature_type,
        cdsStart ? (loc.start - 1) * 3 + cdsStart : loc.start,
        cdsStart ? loc.end * 3 + cdsStart -1 : loc.end,
        score ?: ".",
        strand ?: ".",
        parentId ? "0" : ".",
        attributes.join(";")
    ].join("\t")
}
