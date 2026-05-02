package uk.ac.ebi.interpro

import java.nio.file.Path
import uk.ac.ebi.interpro.Location
import uk.ac.ebi.interpro.LocationFragment
import uk.ac.ebi.interpro.Match
import uk.ac.ebi.interpro.Signature
import uk.ac.ebi.interpro.SignatureLibraryRelease
import uk.ac.ebi.interpro.Site

class SFLD {
    static Map<String, Map<String, Match>> parseOutput(
        Path filePath, 
        Map<String, Map> hmmerMatches
    ) {
        // Parse the output TSV file from the SFLD postprocess bin
        def matches = [:]
        filePath.withReader{ reader ->
            while (true) {
                String sequenceId = null

                for (String line = reader.readLine(); line != null; line = reader.readLine()) {
                    def matcher = line.trim() =~ ~/^Sequence:\s+(\S+)$/
                    if (matcher.find()) {
                        sequenceId = matcher.group(1).trim()
                        break
                    }
                }

                if (sequenceId == null) {
                    break
                }

                assert hmmerMatches.containsKey(sequenceId)
                def seqHmmerMatches = hmmerMatches[sequenceId]
                matches[sequenceId] = parseBlock(reader, sequenceId, seqHmmerMatches)
            }
        }

        return matches
    }

    static Map<String, Match> parseBlock(
        Reader reader,
        String sequenceId,
        Map<String, Map> seqHmmerMatches
    ) {
        SignatureLibraryRelease library = new SignatureLibraryRelease("SFLD", null)
        boolean inDomains = false
        def domains = [:]
        while (true) {
            String line = reader.readLine()?.trim()
            if (!line) break
            if (line == "Domains:") {
                inDomains = true
            } else if (line == "Sites:") {
                inDomains = false
            } else if (line == "//") {
                break
            } else if (inDomains) {
                def fields = line.split(/\t/)
                assert fields.length == 15
                String modelAccession = fields[0]
                Integer hmmStart = fields[4] as Integer
                Integer hmmEnd = fields[5] as Integer
                int aliStart = fields[7] as int
                int aliEnd = fields[8] as int
                Integer envStart = fields[9] as Integer
                Integer envEnd = fields[10] as Integer
                assert seqHmmerMatches.containsKey(modelAccession)
                def hmmerMatch = seqHmmerMatches[modelAccession]
                Location loc = hmmerMatch.locations.find { it.start == aliStart 
                                                            && it.end == aliEnd 
                                                            && it.hmmStart == hmmStart 
                                                            && it.hmmEnd == hmmEnd 
                                                            && it.envelopeStart == envStart 
                                                            && it.envelopeEnd == envEnd }
                assert loc != null
                Match match = domains.computeIfAbsent(modelAccession, k -> {
                    Signature signature = new Signature(modelAccession, library)
                    return new Match(modelAccession, hmmerMatch.evalue, hmmerMatch.score, hmmerMatch.bias, signature)
                });
                match.addLocation(loc)
            } else {
                def fields = line.split(/\s/, 3)
                assert fields.length == 2 || fields.length == 3
                String modelAccession = fields[0]
                String residues = fields[1]
                String description = fields.length == 3 ? fields[2] : null

                Match match = domains.get(modelAccession)
                if (match != null) {
                    Site site = new Site(description, residues)
                    match.addSite(site)                        
                }
            }
        }

        return domains
    }

    static Map<String, Set<String>> getHierarchies(Path filePath) {
        def hierarchies = [:].withDefault { [] as Set }
        filePath.eachLine { line ->
            def nodes = line.split(/\t/).toList()
            nodes.eachWithIndex { node, idx ->
                if (idx > 0) {
                    def ancestors = nodes.subList(0, idx) as Set
                    hierarchies[node].addAll(ancestors)
                }
            }
        }
        return hierarchies
    }
}
