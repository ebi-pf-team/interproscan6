class FastaFile {
    // DNA, RNA, gaps
    static final String NUCLEIC_ALPHABET = "ACGTUN\\.\\s-"
    // 20 standard AAs, Sec, Pyl, any/unknown, Asx, Glx, Xle, gaps
    static final String PROTEIN_ALPHABET = "ACDEFGHIKLMNPQRSTVWYUOXBZJ\\s"

    static Map<String, String> parse(String fastaFilePath) {
        def sequences = [:]
        def md5 = null
        StringBuilder currentSequence = null

        new File(fastaFilePath).eachLine { line ->
            if (line.startsWith(">")) {
                if (currentSequence) {
                    sequences[md5] = currentSequence.toString()
                }
                md5 = line.substring(1).trim()
                currentSequence = new StringBuilder()
            } else {
                currentSequence.append(line.replaceAll("\\s+", ""))
            }
        }
        if (currentSequence) {
            sequences[md5] = currentSequence.toString()
        }
        return sequences
    }
}
