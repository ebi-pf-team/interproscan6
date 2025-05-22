class PirsfDatEntry {
    String modelAccession
    String name
    BigDecimal meanL
    BigDecimal stdL
    BigDecimal minS
    BigDecimal meanS
    BigDecimal stdS
    boolean blast
    Set<String> children

    PirsfDatEntry(String modelAccession) {
        this.modelAccession = modelAccession
    }

    void addName(String name) {
        this.name = name
    }

    void addMeans(String meanL, String stdL, String minS, String meanS, String stdS) {
        this.meanL = new BigDecimal(meanL)
        this.stdL  = new BigDecimal(stdL)
        this.minS  = new BigDecimal(minS)
        this.meanS = new BigDecimal(meanS)
        this.stdS  = new BigDecimal(stdS)
    }

    void addBlastBool(String blast) {
        this.blast = blast != "No"
    }

    void addChildren(String children) {
        this.children = children.split("\\s+")
    }

    static parsePirsfDatFile (String pirsfDatPath) {
        File pirsfDatFile = new File(pirsfDatPath)
        if (!pirsfDatFile.exists()) {
            System.exit(1)
        }

        Map<String, PirsfDatEntry> datEntries = [:]  // modelAccession: entry
        Map<String, String> datChildren = [:]   // subfam: parent

        pirsfDatFile.withReader { reader ->
            String line
            PirsfDatEntry entry = null
            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {
                    if (entry != null) {
                        datEntries.put(entry.modelAccession, entry)
                    }
                    entry = new PirsfDatEntry(line.split("\\s")[0].replace(">","").trim())
                    def childLineMatcher = line =~ ~/^>PIRSF\d+\schild:\s(.*)$/
                    if (childLineMatcher.find()) {
                        entry.addChildren(childLineMatcher.group(1).trim())
                        def children = childLineMatcher.group(1).trim().split("\\s+")
                        for (String subfamily : children) {
                            datChildren.put(subfamily, entry.modelAccession)
                        }
                    }
                }
                else if (line.startsWith("BLAST:")) {
                    entry.addBlastBool(line.replace("BLAST:", "").trim())
                }
                else if (line.matches("\\d+\\.?\\d*\\s+\\d+\\.?\\d*\\s+\\d+\\.?\\d*\\s+\\d+\\.?\\d*\\s+\\d+\\.?\\d*")) {
                    def numbers = line.split("\\s+")
                    entry.addMeans(numbers[0], numbers[1], numbers[2], numbers[3], numbers[4])
                }
                else if (!line.trim().isEmpty()) {
                    entry.addName(line.trim())
                }
            }
            // Add the last entry if it exists
            if (entry != null) {
                datEntries.put(entry.modelAccession, entry)
            }
        }

        return [datEntries, datChildren]
    }
}