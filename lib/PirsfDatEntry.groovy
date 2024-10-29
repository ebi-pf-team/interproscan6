class PirsfDatEntry {
    String modelAccession
    String name
    Float meanL
    Float stdL
    Float minS
    Float meanS
    Float stdS
    boolean blast
    Set<String> children = new HashSet<>()

    PirsfDatEntry(String modelAccession) {
        this.modelAccession = modelAccession
    }

    void addName(String name) {
        this.name = name
    }

    void addMeans(String meanL, String stdL, String minS, String meanS, String stdS) {
        this.meanL = meanL.toFloat()
        this.stdL = stdL.toFloat()
        this.minS = minS.toFloat()
        this.meanS = meanS.toFloat()
        this.stdS = stdS.toFloat()
    }

    void addBlastBool(String blast) {
        this.blast = blast != "No"
    }

    void addChildren(String children) {
        this.children = new HashSet<>(children.split("\\s+").toList())
    }

    static parsePirsfDatFile (String pirsfDatPath) {
        File pirsfDatFile = new File(pirsfDatPath)
        if (!pirsfDatFile.exists()) {
            System.out.println("Could not find PIRSF dat file at ${pirsfDatPath}")
            System.exit(1)  // Added parentheses
        }

        Map<String, PirsfDatEntry> datEntries = new LinkedHashMap<>()  // model id: datEntry
        Map<String, String> datChildren = new LinkedHashMap<>()        // Sub family: Parent family

        pirsfDatFile.withReader { reader ->
            String line
            PirsfDatEntry entry
            while ((line = reader.readLine()) != null) {  // Fixed parentheses
                if (line.startsWith(">")) {
                    if (entry != null) {
                        datEntries.put(entry.modelAccession, entry)
                    }
                    entry = new PirsfDatEntry(line.split("\\s")[0].replace(">","").trim())  // Removed redundant declaration
                    def childLineMatcher = line =~ ~/^>PIRSF\d+\schild:\s(.*)$/
                    if (childLineMatcher.find()) {
                        entry.addChildren(childLineMatcher.group(1).trim())
                        def children = childLineMatcher.group(1).trim().split("\\s+")  // Added 'def' declaration
                        for (String subfamily : children) {
                            datChildren.put(subfamily, entry.modelAccession)  // Fixed map name
                        }
                    }
                }
                else if (line.startsWith("BLAST:")) {  // Fixed parentheses
                    entry.addBlastBool(line.replace("BLAST:", "").trim())
                }
                else {
                    def numLineMatcher = line =~ ~/^(\d+\.?\d+)\s+(\d+\.?\d+)\s+(\d+\.?\d+)\s+(\d+\.?\d+)\s+(\d+\.?\d+)\s*?$/
                    if (numLineMatcher.find()) {  // Fixed parentheses
                        entry.addMeans(
                                numLineMatcher.group(1),
                                numLineMatcher.group(2),
                                numLineMatcher.group(3),
                                numLineMatcher.group(4),
                                numLineMatcher.group(5)
                        )
                    }
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