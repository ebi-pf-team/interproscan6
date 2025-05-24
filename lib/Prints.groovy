class Prints {
    String modelName
    String modelId
    Double evalue
    String graphscan
    Integer locationStart
    Integer locationEnd
    Double pvalue
    Double score
    Integer motifNumber
    Integer motifCount

    Prints(String modelName, String modelId, Double evalue) {
        this.modelName = modelName
        this.modelId = modelId
        this.evalue = evalue
    }

    static Prints buildLocationMatch(
            Prints printsMatch,
            int start,
            int end,
            double pvalue,
            double score,
            int motifNum,
            int motifCount
    ) {
        Prints newMatch = new Prints(printsMatch.modelName, printsMatch.modelId, printsMatch.evalue)
        newMatch.graphscan = printsMatch.graphscan
        newMatch.locationStart = start
        newMatch.locationEnd = end
        newMatch.pvalue = pvalue
        newMatch.score = score
        newMatch.motifNumber = motifNum
        newMatch.motifCount = motifCount
        return newMatch
    }
}

class HierarchyEntry { // represents a row in the HierarchyDB
    String modelId
    String modelAccession
    Double evalueCutoff
    int minMotifCount
    boolean isDomain = false
    String[] siblingsIds = []

    HierarchyEntry(String modelId, String modelAccession, Double evalueCutoff, int minMotifCount) {
        this.modelId = modelId
        this.modelAccession = modelAccession
        this.evalueCutoff = evalueCutoff
        this.minMotifCount = minMotifCount
    }

    static parseHierarchyDbFile(String hierarchyDb) {  // parser the hierarchyDb in to a map of ModelID: hierarchy entry
        Map<String, HierarchyEntry> hierarchyMap = new LinkedHashMap<>()
        File hierarchyFile = new File(hierarchyDb)
        if (!hierarchyFile.exists()){
            System.out.println("Could not find Hierarchy DB for PRINTS")
            System.exit 1
        }
        hierarchyFile.withReader { reader ->
            String line
            while ((line = reader.readLine()) != null) {
                if (!line.startsWith("#")) {
                    def row = line.trim().split("\\|")
                    if (row.length >= 3) {
                        // e.g. G6PDHDRGNASE|PR00079|1e-04|0|
                        final String modelId = row[0].trim()
                        final String modelAccession = row[1].trim()
                        final Double evalueCutoff = Double.parseDouble(row[2].trim())
                        final int minMotifCount = Integer.parseInt(row[3].trim())
                        HierarchyEntry hierarchyEntry = new HierarchyEntry(modelId, modelAccession, evalueCutoff, minMotifCount)
                        if (row.length > 4) {
                            // e.g. DHBDHDRGNASE|PR01397|1e-04|0|SDRFAMILY or TYROSINASE|PR00092|1e-04|0|*
                            String siblingsOrDomain = row[4].trim()
                            if (siblingsOrDomain == "*") {
                                hierarchyEntry.setDomainStatus()
                            } else if (siblingsOrDomain.size() > 0) {
                                hierarchyEntry.addSiblings(siblingsOrDomain)
                            }
                        }
                        hierarchyMap.put(modelId, hierarchyEntry)
                    }
                }
            }
        }
        return hierarchyMap
    }

    void setDomainStatus() {
        this.isDomain = true
    }

    void addSiblings(String siblingsString) {
        String[] siblingsIds = siblingsString.split("\\,")
        this.siblingsIds = siblingsIds
    }
}
