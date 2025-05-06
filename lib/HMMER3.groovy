class HMMER3 {
    static parseOutput(String filePath, String memberDb) {
        File file = new File(filePath)
        String line
        String queryName
        Integer queryLength
        String queryAccession
        String targetId
        def hits = [:].withDefault { [:] }
        file.withReader { reader ->
            while (true) {
                // Move to the next Query block
                while (queryName == null) {
                    line = reader.readLine()
                    if (line == null) {
                        break
                    }
                    def matcher = line =~ ~/^Query:\s*(.+)\s+\[\w=(\d+)\]$/
                    if (matcher.find()) {
                        queryName = matcher.group(1).trim()
                        queryLength = matcher.group(2).toInteger()
                        break
                    }
                }

                if (queryName == null) {
                    // We reached the end of the file
                    break
                }

                // Move to the sequence top hits list and skip the header
                while (!line.startsWith("Scores for complete sequences")) {
                    def matcher = line =~ ~/^Accession:\s*(.+)$/
                    if (matcher.find()) {
                        queryAccession = matcher.group(1).trim()
                    }

                    line = reader.readLine()
                }
                line = reader.readLine()
                line = reader.readLine()
                line = reader.readLine()

                // Fallback to query name if HMM doesn't have an ACC field
                queryAccession ?= queryName

                // Parse the sequence top hits
                boolean isIncluded = true
                while ((true)) {
                    line = reader.readLine().trim()
                    if (line.isEmpty() || 
                        line.contains("[No hits detected that satisfy reporting thresholds]")) {
                        break
                    } else if (line.startsWith("------ inclusion threshold")) {
                        isIncluded = false
                    }

                    def fields = line.split(/\s+/)
                    Match match = new Match(
                        queryAccession,
                        Double.parseDouble(fields[0]),
                        Double.parseDouble(fields[1]),
                        Double.parseDouble(fields[2]),
                        new Signature(queryAccession)
                    )
                    match.signature.signatureLibraryRelease.library = memberDb
                    match.included = isIncluded
                    targetId = fields[8].trim()
                    hits[targetId][queryAccession] = match
                }

                // Move to the hits section
                while (!line.startsWith("Domain annotation for")) {
                    line = reader.readLine().trim()
                }

                // Parse hits
                while (true) {
                    // Find next target
                    targetId = null
                    while (true) {
                        if (line.startsWith(">>")) {
                            def matcher = line =~ ~/>>\s(\S+)/
                            assert matcher.find() == true
                            targetId = matcher.group(1)
                            break
                        } else if (line == "//") {
                            break
                        }
                        line = reader.readLine().trim()
                    }

                    if (targetId == null) {
                        // Reached the end of the query block
                        break
                    }

                    // At this point we need to check whether we have individual domains
                    line = reader.readLine()
                    if (line.trim().startsWith("[No individual domains")) {
                        // We do not
                        continue
                    }

                    // Move the first domain
                    line = reader.readLine()

                    // Parse domain hits
                    while (!(line = reader.readLine().trim()).isEmpty()) {
                        def fields = line.split(/\s+/)
                        assert fields.size() == 16
                        Location location = new Location(
                            Integer.parseInt(fields[9]),
                            Integer.parseInt(fields[10]),
                            fields[6].toInteger(),
                            fields[7].toInteger(),
                            queryLength,
                            fields[8],
                            fields[12].toInteger(),
                            fields[13].toInteger(),
                            Double.parseDouble(fields[5]),
                            Double.parseDouble(fields[2]),
                            Double.parseDouble(fields[3])
                        )
                        location.included = fields[1] == "!"
                        hits[targetId][queryAccession]?.addLocation(location)
                    }

                    // Move to domain alignments
                    while (true) {
                        line = reader.readLine().trim()
                        if (line == "Alignments for each domain:") {
                            break
                        }
                    }

                    // Find the alignment for each domain
                    int domainIndex = 0
                    String queryAlignment = ""
                    String targetAlignment = ""
                    while (true) {
                        line = reader.readLine().trim()
                        if (line.isEmpty()) {
                            continue
                        } else if (line =~ /==\s+domain\s\d+/) {
                            // New domain
                            if (domainIndex > 0) {
                                assert queryAlignment.length() > 0
                                assert targetAlignment.length() > 0
                                hits[targetId][queryAccession].setAlignments(
                                    domainIndex - 1,
                                    queryAlignment,
                                    targetAlignment
                                )
                            }

                            domainIndex++
                            queryAlignment = ""
                            targetAlignment = ""
                        } else if (line.startsWith(">>") || 
                                   line == "Internal pipeline statistics summary:") {
                            // Either new target or end of results for current query

                            assert queryAlignment.length() > 0
                            assert targetAlignment.length() > 0
                            hits[targetId][queryAccession].setAlignments(
                                domainIndex - 1,
                                queryAlignment,
                                targetAlignment
                            )

                            break
                        } else {
                            // Read next block of alignments
                            def blocks = []
                            while (!line.isEmpty() || blocks.size() < 4) {
                                blocks.add(line)
                                line = reader.readLine().trim()
                            }
                            /* 
                                Keep last four lines as the domain block may have five lines 
                                if the profile has a consensus structure 
                                or reference line annotation (reported first)
                            */
                            blocks = blocks[-4..-1]
                            // Consensus of query profile
                            def fields = blocks[0].split(/\s+/)
                            queryAlignment += fields[2]
                            /*
                                Target sequence, with dashes (-) indicate deletions 
                                in the target sequence with respect to the profile
                            */
                            fields = blocks[2].split(/\s+/)
                            targetAlignment += fields[2]
                        }
                    }

                }

                queryName = null;
                queryAccession = null
            }                
        }

        // Remove matches without locations, 
        hits.entrySet().removeIf { Map.Entry<String,Map<String,Match>> targetEntry ->
            targetEntry.value.entrySet().removeIf { Map.Entry<String,Match> matchEntry ->
                matchEntry.value.locations.isEmpty()
            }
            targetEntry.value.isEmpty()
        }

        return hits
    }

    static splitByLocation(hmmerMatches) {
        // This method is used by CATH-Gene3D and CATH-FunFam
        def results = [:].withDefault { [:] }
        hmmerMatches.each { sequenceId, matches ->
            matches.values().each { m1 ->
                m1.locations.each { loc ->
                    Match m2 = new Match(
                        m1.modelAccession, 
                        m1.evalue,
                        m1.score, 
                        m1.bias
                    )
                    m2.addLocation(loc)
                    String key = "${m2.modelAccession}-${loc.envelopeStart}-${loc.envelopeEnd}"
                    results[sequenceId][key] = m2
                }
            }
        }
        return results
    }
}
