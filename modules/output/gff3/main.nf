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

    gff3File.withWriter { gff3Writer ->
        gff3Writer.writeLine("##gff-version 3.1.26")
        gff3Writer.writeLine("##interproscan-version ${interproscan_version}")

        def tempFastaFile = new File("temp.fasta")
        tempFastaFile.withWriter { fastaWriter ->
            fastaWriter.writeLine("##FASTA")

            Set<String> seenNucleicMd5s = new HashSet<>()

            Pattern esl_pattern = Pattern.compile(/^source=[^"]+\s+coords=(\d+)\.\.(\d+)\s+length=\d+\s+frame=(\d+)\s+desc=.*$/)

            matches_files.each { matchFile ->
                matchFile = new File(matchFile.toString())
                Map proteins = new ObjectMapper().readValue(matchFile, Map)

                if (nucleic) {
                    // Collect all protein MD5s
                    def allProteinMd5s = proteins.keySet() as List

                    // Retrieve nucleic-protein associations in batches
                    def nucleicToProteinMd5 = [:].withDefault { [] as Set }
                    allProteinMd5s.collate(1000).each { batch ->
                        def partial = db.groupProteins(batch)
                        partial.each { ntMd5, protMd5s ->
                            nucleicToProteinMd5[ntMd5].addAll(protMd5s)
                        }
                    }

                    // Retrieve nucleotide sequences in batches
                    def allNtMd5s = nucleicToProteinMd5.keySet() as List
                    def ntSeqDataMap = [:].withDefault { [] }
                    allNtMd5s.collate(1000).each { batch ->
                        def partial = db.nucleicMd5ToNucleicSeqs(batch)
                        ntSeqDataMap.putAll(partial)
                    }

                    // Build all (protein, nucleotide) pairs
                    def allPairs = []
                    def seenPairKeys = new HashSet<String>()
                    nucleicToProteinMd5.each { ntMd5, protMd5s ->
                        protMd5s.each { pMd5 ->
                            def key = "${pMd5}\t${ntMd5}"
                            if (!seenPairKeys.contains(key)) {
                                seenPairKeys.add(key)
                                allPairs << [pMd5, ntMd5]
                            }
                        }
                    }

                    // Retrieve ORF data in batches
                    def orfDataMap = [:].withDefault { [:].withDefault { [] } }
                    allPairs.collate(500).each { batch ->
                        def partial = db.getOrfSeqs(batch)
                        partial.each { ntMd5, protMap ->
                            protMap.each { protMd5, dataList ->
                                orfDataMap[ntMd5][protMd5].addAll(dataList)
                            }
                        }
                    }

                    nucleicToProteinMd5.each { String nucleicMd5, Set<String> proteinMd5s ->
                        if (seenNucleicMd5s.contains(nucleicMd5)) {
                            return
                        }
                        seenNucleicMd5s.add(nucleicMd5)

                        def seqData = ntSeqDataMap[nucleicMd5]
                        seqData.each { seq ->
                            String seqId = seq.id
                            String sequence = seq.sequence.trim()
                            int seqLength = sequence.length()
                            gff3Writer.writeLine("##sequence-region ${seqId} 1 ${seqLength}")

                            proteinMd5s.each { String proteinMd5 ->
                                def proteinSeqData = orfDataMap[nucleicMd5][proteinMd5]
                                proteinSeqData.each { row ->
                                    def matcher = esl_pattern.matcher(row.description)
                                    assert matcher.matches()
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

                                    proteins[proteinMd5].each { modelAcc, match->
                                        match = Match.fromMap(match)
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
                } else {
                    def allMd5s = proteins.keySet() as List
                    def md5ToSeqData = [:]
                    allMd5s.collate(1000).each { batch ->
                        def result = db.proteinMd5ToProteinSeqs(batch)
                        md5ToSeqData.putAll(result)
                    }

                    proteins.each { String proteinMd5, Map matchesMap ->
                        def seqData = md5ToSeqData[proteinMd5]
                        String sequence = seqData[0].sequence.trim()
                        int seqLength = sequence.length()

                        seqData.each { row ->
                            gff3Writer.writeLine("##sequence-region ${row.id} 1 ${seqLength}")

                            matchesMap.each { modelAcc, match ->
                                match = Match.fromMap(match)
                                
                                match.locations.each { Location loc ->
                                    gff3Writer.writeLine(proteinFormatLine(row.id, match, loc, null, null, null))
                                }
                            }

                            fastaWriter.writeLine(">${row.id}")
                            fastaWriter.writeLine("${sequence.replaceAll(/(.{60})/, '$1\n')}")   
                        }         
                    }
                }
            }
        }

        tempFastaFile.withReader { fastaReader ->
            fastaReader.eachLine { line ->
                gff3Writer.writeLine(line)
            }
        }
        tempFastaFile.delete()
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
        case "TMHMM":
        case "DeepTMHMM":
            feature_type = match.signature.accession.toUpperCase() == "TRANSMEMBRANE ALPHA HELIX" ? "transmembrane_helix" :
                match.signature.accession.toUpperCase() == "TRANSMEMBRANE BETA BARREL" ? "transmembrane_polypeptide_region" :
                match.signature.accession.toUpperCase() == "PERIPLASMIC DOMAIN" ? "non_cytoplasmic_polypeptide_region" : "signal_peptide"
            break
        case "Phobius":
            feature_type = match.signature.description.toUpperCase() == "CYTOPLASMIC DOMAIN" ? "cytoplasmic_polypeptide_region" :
                match.signature.description.toUpperCase() == "NON CYTOPLASMIC DOMAIN" ? "non_cytoplasmic_polypeptide_region" :
                match.signature.description.toUpperCase() == "TRANSMEMBRANE REGION" ? "transmembrane_polypeptide_region" : "signal_peptide"
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
