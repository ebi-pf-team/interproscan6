package uk.ac.ebi.interpro
import java.util.regex.Pattern
import com.fasterxml.jackson.core.JsonFactory
import com.fasterxml.jackson.databind.ObjectMapper
import uk.ac.ebi.interpro.Location
import uk.ac.ebi.interpro.Match
import uk.ac.ebi.interpro.SeqDB


class ProcessOutputGFF3 {
    static void run(List<String> inputPaths, String databasePath, boolean isNucleic, String iprscanVersion, String outputPath) {
        SeqDB db = new SeqDB(databasePath)

        new File(outputPath).withWriter { gff3Writer ->
            gff3Writer.writeLine("##gff-version 3.1.26")
            gff3Writer.writeLine("##interproscan-version ${iprscanVersion}")

            def tempFastaFile = new File("${outputPath}.fasta")
            tempFastaFile.withWriter { fastaWriter ->
                fastaWriter.writeLine("##FASTA")

                Set<String> seenNucleicMd5s = new HashSet<>()

                Pattern esl_pattern = Pattern.compile(/^source=[^"]+\s+coords=(\d+)\.\.(\d+)\s+length=\d+\s+frame=(\d+)\s+desc=.*$/)

                ObjectMapper mapper = new ObjectMapper()
                inputPaths.each { matchFile ->
                    matchFile = new File(matchFile.toString())
                    Map proteins = mapper.readValue(matchFile, Map)

                    if (isNucleic) {
                        def (nucleicToProteinMd5, ntSeqDataMap, orfDataMap) = 
                            db.retrieveAllNucleicSequenceData(proteins.keySet() as List)

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

    static String proteinFormatLine(String seqId, Match match, Location loc, 
                                    String parentId, Integer cdsStart, String strand) {
        String memberDb = match.signature.signatureLibraryRelease.library

        def goTerms = []
        if (memberDb == "PANTHER" && match.treegrafter?.goXRefs){
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
                feature_type = ["DOMAIN", "REPEAT"].contains(match.signature.type?.toUpperCase()) ? "polypeptide_domain" : "polypeptide_region"
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
        if (match.source == "InterPro-N") {
            score = loc.score
        } else if (["CDD", "PRINTS"].contains(memberDb)) {
            score = match.evalue
        } else if (["HAMAP", "PROSITE profiles"].contains(memberDb)) {
            score = loc.score
        } else if (["SignalP-Prok", "SignalP-Euk"].contains(memberDb)) {
            score = loc.pvalue
        } else {
            score = loc.evalue
        }

        def name
        def alias = null
        if (match.signature.name) {
            name = match.signature.name
            alias = match.signature.accession
        } else {
            name = match.signature.accession
        }

        String interproAccession = match.signature.entry?.accession
        def attributes = [
            "Name=${name}",
            alias ? "Alias=${alias}" : null,
            parentId ? "Parent=${parentId}" : null,
            interproAccession ? "Dbxref=InterPro:${interproAccession}" : null,
            (match.source == "InterPro-N") ? "Dbxref=${memberDb}:${match.signature.accession}" : null,
            goTerms ? "Ontology_term=${goTerms.collect { it.id }.join(',')}" : null,
            match.signature.type ? "type=${match.signature.type}" : null,
            "representative=${loc.representative}",
        ].findAll { it }

        return [
            seqId,
            match.source,
            feature_type,
            cdsStart ? (loc.start - 1) * 3 + cdsStart : loc.start,
            cdsStart ? loc.end * 3 + cdsStart -1 : loc.end,
            score ?: ".",
            strand ?: ".",
            parentId ? "0" : ".",
            attributes.join(";")
        ].join("\t")
    }
}
