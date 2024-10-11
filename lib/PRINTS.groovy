class PRINTS {
    static parseOutput(
            String printsOutput,
            String hierarchyDb
    ) {
        // Build up a map of the fingerprint hierarchies
        Map<String, Map<String, Object>> hierarchyMap = new LinkedHashMap<>()
        File hierarchyFile = new File(hierarchyDb)
        hierarchyFile.withReader { reader ->
            String line
            while ((line = reader.readLine()) != null) {
                if (!line.startsWithAny("#", "/")) {
                    def row = line.split("\\|")
                    if (row.length >= 3) {
                        String motifName = row[0]
                        String motifAccession = row[1]
                        Double evalueCutoff = Double.parseDouble(row[2])
                        int minMotifCount = Integer.parseInt(row[3])
                        hierarchyMap.put(
                                motifName,
                                [
                                        "motifAccession": motifAccession,
                                        "cutOff": evalueCutoff,
                                        "minMotifCount": minMotifCount
                                ]
                        )
                    }
                }
            }
        }

        // parse the Prints Output into Match instances
        Map<String, Map<String, Match>> hits = new LinkedHashMap<>()
        Map<String, String> name2accession = new LinkedHashMap<>()
        File printsFile = new File(printsOutput)
        printsFile.withReader { reader ->
            String line
            String queryAccession = ""
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("Sn;")) {  // Get the query sequence id
                    queryAccession = line.replaceFirst("Sn; ", "").trim()
                    name2accession = new LinkedHashMap<>()
                }
                // Do not retrieve the motif Description from the 1TBH line - this is retrieved in the XREFS subworkflow
                // And retrieve the motifAcc from the heirarchyDB instead of the 1TBH line
                else if (line.startsWith("2TBH")) {  // Used to create Match instances
                    // Line: 2TBH  MotifName  NumMotifs  SumId  AveId  ProfScore  Ppvalue  Evalue  GraphScan
                    def matcher = line =~ ~/^2TBH\s+(.+?)\s+(\d)\s+of\s+\d\s+\d+\.?\d*\s+\d+\.?\d*\s+\d+\s+.*?\s+(.*?)\s+([\.iI]+)\s+$/
                    if (matcher.find()) {
                        String motifName = matcher.group(1).trim()
                        int numOfMotifs = Integer.parseInt(matcher.group(2).trim())
                        Double evalue = Double.parseDouble(matcher.group(3).trim())
                        String graphScan = matcher.group(4).trim()

                        if (hierarchyMap.containsKey(motifName)) {
                            String motifAccession = hierarchyMap[motifName]["motifAccession"]
                            int minMotifCount = hierarchyMap[motifName]["minMotifCount"] as int
                            Double cutoff = hierarchyMap[motifName]["cutOff"] as Double

                            if (evalue <= cutoff && numOfMotifs > minMotifCount) {
                                if (!hits.containsKey(queryAccession)) {
                                    hits.put(queryAccession, new LinkedHashMap<>())
                                }
                                hits[queryAccession].put(motifAccession, new Match(motifAccession, evalue, graphScan))
                                name2accession.put(motifName, motifAccession)
                            }
                        }
                    }
                }
                else if (line.startsWith("3TBH")) {  // Used to create Location instances
                    // Line: 3TBH MotifName NoOfMotifs IdScore PfScore Pvalue Sequence Len Low Pos High
                    def matcher = line =~ ~/^3TBH\s+(.+?)\s+(\d+)\s+of\s+\d+\s+([\d\.]+)\s+\d+\s+(.+?)\s+([\w#]+)\s+(\d+)\s+\d+\s+(-?\d+)\s*\d\s*$/
                    if (matcher.find()) {
                        // Check motif passed the e-value cutoff when parsing line 2TBH - else SKIP!
                        String motifName = matcher.group(1).trim()
                        String motifAccession = name2accession[motifName]
                        if (hits.containsKey(queryAccession)) {
                            if (hits[queryAccession].containsKey(motifAccession)) {
                                int motifNumber = Integer.parseInt(matcher.group(2).trim())
                                Float idScore = Float.parseFloat(matcher.group(3).trim())
                                Double pvalue = Double.parseDouble(matcher.group(4).trim())
                                String sequence = matcher.group(5).trim()
                                int motifLength = Integer.parseInt(matcher.group(6).trim())
                                int position = Integer.parseInt(matcher.group(7).trim())

                                int end = position + motifLength - 1
                                if (position < 1) {
                                    position = 1
                                }

                                int indexCheck = 0
                                if (sequence.endsWith("#")) {
                                    motifLength = sequence.length()
                                    indexCheck = motifLength - 1
                                    while (sequence[indexCheck] == "#") {
                                        indexCheck -= 1
                                    }
                                }
                                end = end - (motifLength - indexCheck) + 1

                                Location location = new Location(position, end, pvalue, idScore, motifNumber)
                                hits[queryAccession][motifAccession].addLocation(location)
                            }
                        }
                    }
                }
            }
        }

        return hits
    }
}
