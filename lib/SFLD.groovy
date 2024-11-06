class SFLD {
    static Map<String, Set<String>> parseHierarchy(String filePath) {
        def hierarchy = [:]
        File file = new File(filePath)
        file.withReader{ reader ->
            for (String line = reader.readLine(); line != null; line = reader.readLine()) {
                def fields = line.split(":")
                assert fields.size() == 2
                String childAccession = fields[0]
                def parents = fields[1].trim().split(/\s+/)
                hierarchy[childAccession] = (parents*.trim()) as Set
            }
        }

        return hierarchy
    }

    static Map<String, Map<String, Match>> parseOutput(String outputFilePath) {
        def matches = [:]
        File file = new File(outputFilePath)
        file.withReader{ reader ->
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

                matches[sequenceId] = SFLD.parseBlock(reader)
            }
        }

        return matches
    }

    static Map<String, Match> parseBlock(Reader reader) {
        boolean inDomains = false
        def domains = [:]
        while (true) {
            String line = reader.readLine().trim()
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
                Double seqEvalue = fields[1] as Double
                Double seqScore = fields[2] as Double
                Double seqBias = fields[3] as Double
                Integer hmmStart = fields[4] as Integer
                Integer hmmEnd = fields[5] as Integer
                Double domScore = fields[6] as Double
                int aliStart = fields[7] as int
                int aliEnd = fields[8] as int
                Integer envStart = fields[9] as Integer
                Integer envEnd = fields[10] as Integer
                Double domCEvalue = fields[11] as Double
                Double domIEvalue = fields[12] as Double
                Double accuracy = fields[13] as Double
                Double domBias = fields[14] as Double

                if (aliStart <= aliEnd && hmmStart <= hmmEnd && envStart <= envEnd) {
                    Match match = domains.computeIfAbsent(modelAccession, k -> {
                        return new Match(modelAccession, seqEvalue, seqScore, seqBias)
                    });
                    Location location = new Location(aliStart, aliEnd, hmmStart, hmmEnd, null, null,
                            envStart, envEnd, domIEvalue, domScore, domBias)
                    match.addLocation(location)
                }
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
}
