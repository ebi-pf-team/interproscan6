class PIRSFDatEntry {
    String modelAccession
    String name
    Float mean_l
    Float std_l
    Float min_s
    Float mean_s
    Float std_s
    boolean blast
    Set<String> childern = []

    PIRSFDatEntry(String modelAccession) {
        this.modelAccession = modelAccession
    }

    void addName(String name) {
        this.name = name
    }

    void addMeans(Float mean_l, Float std_l, Float min_s, Float mean_s, Float std_s) {
        this.mean_l = mean_l
        this.std_l = std_l
        this.min_s = min_s
        this.mean_s = mean_s
        this.std_s = std_s
    }

    void addBlastBool(String blast) {
        this.blast = blast != "No"
    }

    void addChildren(String children) {
        this.childern = children.split("\\s+")
    }
}
