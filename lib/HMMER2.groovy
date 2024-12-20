class HMMER2 {
    static parseOutput(String filePath, Map<String, Integer> hmmLengths) {
        File file = new File(filePath)
        String line
        String querySequence
        def hits = [:].withDefault { [:] }
        SignatureLibraryRelease library = new SignatureLibraryRelease("SUPERFAMILY", null)

        file.withReader { reader ->
            while (true) {
                // Move to the next Query block
                while (querySequence == null) {
                    line = reader.readLine()
                    if (line == null) {
                        break
                    }
                    def matcher = line =~ ~/^Query sequence:\s*(.+)$/
                    if (matcher.find()) {
                        querySequence = matcher.group(1).trim()
                        break
                    }
                }

                if (querySequence == null) {
                    // We reached the end of the file
                    break
                }

                // Move the beginning of the table of sequence hits
                while (!line.startsWith("Scores for sequence family")) {
                    line = reader.readLine()
                }

                def sequenceHits = [:]
                def domainsPerHit = [:]
                while (true) {
                    line = reader.readLine().trim()

                    if (line.startsWith("Model") || 
                        line.startsWith("--------") || 
                        line.isEmpty()) {
                        // Header of empty line
                        continue
                    } else if (line.contains("[no hits above thresholds]")) {
                        // No hits found
                        break
                    } else if (line.startsWith("Parsed for domains")) {
                        // Beginning of the domain hits
                        break
                    } 

                    def fields = line.split(/\s+/)
                    assert fields.size() >= 4
                    /*
                    Model          Description                              Score    E-value  N 
                    --------       -----------                              -----    ------- ---
                    SM00870                                                 462.1   2.8e-134   1

                    ^ If there is a description, it's going to be splitted on spaces as well.
                        We need the first field (accession) and the last three
                    */
                    String modelAccession = fields[0]
                    def tail = fields.takeRight(3)
                    assert tail.size() == 3
                    double score = Double.parseDouble(tail[0])
                    double evalue = Double.parseDouble(tail[1])
                    int numDomains = tail[2].toInteger()

                    Match match = new Match(modelAccession, evalue, score)
                    match.signature = new Signature(modelAccession, library)
                    sequenceHits[modelAccession] = match
                    domainsPerHit[modelAccession] = numDomains
                }

                if (!sequenceHits.isEmpty()) {
                    // We should be at the beginning of the domain hits list
                    assert line.startsWith("Parsed for domains")

                    while (true) {
                        line = reader.readLine().trim()
                        
                        if (line.startsWith("Model") || 
                            line.startsWith("--------") || 
                            line.isEmpty()) {
                            // Header of empty line
                            continue
                        } else if (line == "//") {
                            // End of record
                            break
                        }

                        def fields = line.split(/\s+/)
                        assert fields.size() == 10

                        String modelAccession = fields[0]
                        int start = fields[2].toInteger()
                        int end = fields[3].toInteger()
                        String seqBounds = fields[4]
                        Integer hmmStart = Integer.parseInt(fields[5])
                        Integer hmmEnd = Integer.parseInt(fields[6])
                        String hmmBounds = fields[7]
                        Double score = Double.parseDouble(fields[8])
                        Double evalue = Double.parseDouble(fields[9])
                        Integer hmmLength = hmmLengths[modelAccession]
                        assert hmmLength != null
                        
                        Location loc = new Location(start, end, hmmStart, hmmEnd, hmmLength, hmmBounds,
                                                    null, null, evalue, score, null)
                        
                        Match match = sequenceHits[modelAccession]
                        assert match != null
                        match.addLocation(loc)
                    }
                }

                sequenceHits.each { modelAccession, match -> {
                    // Check that all domain hits have been converted in locations
                    int n = domainsPerHit[modelAccession]
                    assert match.locations.size() == n

                    hits[querySequence][modelAccession] = match
                }}

                querySequence = null
            }                
        }

        return hits
    }

    static parseHMM(String filePath) {
        def hmmLengths = [:]
        String modelAccession = null
        Integer length = null
        new File(filePath).eachLine { line -> 
            if (line.startsWith("ACC ")) {
                def fields = line.split(/\s+/)
                assert fields.size() == 2
                modelAccession = fields[1]
            } else if (line.startsWith("LENG ")) {
                def fields = line.split(/\s+/)
                assert fields.size() == 2
                length = Integer.parseInt(fields[1])
            } else if (line.startsWith("//")) {
                assert modelAccession != null
                assert length != null
                hmmLengths[modelAccession] = length
                modelAccession = null
                length = null
            }
        }

        return hmmLengths
    }
}
