import Match

class HMMER3 {
    static final String queryNamePattern = ~/^Query:\s*(.+)\s+\[\w=(\d+)\]$/
    static final String queryAccPattern = ~/^Accession:\s*(.+)$/
    static final String noHitsPattern = "[No hits detected that satisfy reporting thresholds]"
    static final String noDomainsPattern = "[No individual domains that satisfy reporting thresholds"
    static final String domainsTableHeaderPattern = " ---   ------ ----- --------- ---------"
    static final String targetIdPattern = ~/>>\s(\S+)/

    static void parseOutput(String filePath) {
        File file = new File(filePath)
        String line
        String queryName
        Integer queryLength
        String queryAccession
        String targetId
        def hits = [:]

        file.withReader { reader ->
            while (true) {
                // Move to the next Query block
                while (queryName == null) {
                    line = reader.readLine()
                    if (line == null) {
                        break
                    }
                    def matcher = line =~ this.queryNamePattern
                    if (matcher.find()) {
                        queryName = matcher.group(1)
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
                    def matcher = line =~ this.queryAccPattern
                    if (matcher.find()) {
                        queryAccession = matcher.group(1)
                    }

                    line = reader.readLine()
                }
                line = reader.readLine()
                line = reader.readLine()
                line = reader.readLine()

                // Parse the sequence top hits
                while ((true)) {
                    line = reader.readLine().trim()
                    if (line.isEmpty() || line.contains(this.noHitsPattern)) {
                        break
                    }

                    def fields =line.split(/\s+/)
                    Match match = new Match(
                        queryAccession ?: queryName,
                        Double.parseDouble(fields[0]),
                        Double.parseDouble(fields[1]),
                        Double.parseDouble(fields[2]),
                    )
                    targetId = fields[8]

                    if (!hits.containsKey(targetId)) {
                        hits[targetId] = [:]
                    } 

                    hits[targetId][match.modelAccession] = match
                }

                // Move to the hits section
                while (!line.startsWith("Domain annotation for")) {
                    line = reader.readLine().trim()
                }

                while (true) {
                    line = reader.readLine().trim()
                    if (line.isEmpty() || line.contains(this.noDomainsPattern)) {
                        break
                    }
                }

                // Move to the hits section
                while (line.trim() != "//") {
                    if (line.contains(this.noHitsPattern)) {
                        // No hits found: move to next query
                        while ((line = reader.readLine()) != "//") {
                            continue
                        }
                    } else if (line.startsWith("Domain annotation for")) {
                        // Beginning of the hits section

                        // Find the target
                        while (!(line = reader.readLine()).startsWith(">>")) {
                            continue
                        }

                        def matcher = line =~ this.targetIdPattern
                        assert matcher.find() == true
                        targetId = matcher.group(1)

                        // Move to the beginning of the domain table
                        boolean hasDomains = true
                        while (true) {
                            line = reader.readLine()
                            if (line.startsWith(this.noDomainsPattern)) {
                                hasDomains = false
                                break
                            } else if (line.startsWith(this.domainsTableHeaderPattern)) {
                                break
                            }
                        }

                        if (! hasDomains) {
                            continue
                        }

                        // Parse individual domains
                        while (!(line = reader.readLine().trim()).isEmpty()) {
                            def fields = line.split(/\s+/)
                            assert fields.size() == 16
                            // Match match = new Match(queryAccession ?: queryName)
                            // match.addLocation(
                            //     Integer.parseInt(fields[9]),
                            //     Integer.parseInt(fields[10]),
                            //     fields[6].toInteger(),
                            //     fields[7].toInteger(),
                            //     0,  // TODO: hmmLength
                            //     fields[8],
                            //     fields[12].toInteger(),
                            //     fields[13].toInteger(),
                            //     // evalue
                            //     // score
                            //     // bias
                            //     // alignment
                            //     // feature
                            //     // fragment
                            // )
                        
                        }
                    } else {
                        line = reader.readLine()
                    }
                }

                // We reached the end of the results for the current query
                // println(queryName)

                queryName = null;
                queryAccession = null
            }
            
            // while ((line = reader.readLine()) != null) {
            //     String query = null;
                
                    
                
            
            // }
                
        }
        

        // file.eachLine { line ->
        //     if (line.startsWith(">> ")) {
        //         String targetName = line.substring(3)
        //     } else if (line.startsWith("== domain")) {

        //     }
        // }
    }    
}
