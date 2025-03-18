class FastaFile {
    // DNA, RNA, gaps
    static final String NUCLEIC_ALPHABET = "ACGTUN\\.\\s-"
    // 20 standard AAs, Sec, Pyl, any/unknown, Asx, Glx, Xle, gaps
    static final String PROTEIN_ALPHABET = "ACDEFGHIKLMNPQRSTVWYUOXBZJ\\.\\s-"

    static Map<String, String> parse(String fastaFilePath) {
        def sequences = [:]
        def md5 = null
        def currentSequence = null
        new File(fastaFilePath).eachLine { line ->
            if (line.startsWith(">")) {
                if (currentSequence) {
                    sequences[md5] = currentSequence
                }
                md5 = line.substring(1).trim()
                currentSequence = ""
            } else {
                currentSequence += line.trim()
            }
        }
        if (currentSequence) {
            sequences[md5] = currentSequence
        }
        return sequences
    }

    static ArrayList validate(String fastaFilePath, boolean isNucleic, Map appsConfig, List<String> appsToRun) {
        def cmd = null
        String errorMsg = ""

        // Check for generally illegal chars
        String alphabet = isNucleic ? NUCLEIC_ALPHABET : PROTEIN_ALPHABET
        cmd = "grep -v '^>' $fastaFilePath | grep -En '[^$alphabet${alphabet.toLowerCase()}]'"
        String grepOut = searchFile(cmd)
        if (grepOut) {
            errorMsg += "Forbidden characters found in $fastaFilePath on lines: $grepOut\n"
        }

        // Check for application-specific forbidden chars
        appsToRun.each{ app ->
            def forbiddenChars = appsConfig[app].invalid_chars ?: ""
            forbiddenChars = forbiddenChars.toSet()  // remove duplicates
            forbiddenChars.each { String forbiddenChar ->
                cmd = "grep -En '^[^>].*[$forbiddenChar${forbiddenChar.toLowerCase()}]' $fastaFilePath"
                grepOut = searchFile(cmd)
                if (grepOut) {
                    errorMsg += "$app forbidden character '$forbiddenChar' found in $fastaFilePath on lines: $grepOut\n"
                }
            }
        }

        // Count the number of sequences in the input FASTA file
        cmd = ["grep", "-c", "^>", fastaFilePath]
        def process = cmd.execute()
        process.waitFor()
        Integer numSequences = process.in.text.trim() as Integer

        return [errorMsg, numSequences]
    }

    static String searchFile(def cmd) {
        try {
            def process = ["bash", "-c", cmd].execute()
            process.waitFor()
            def output = process.in.text.trim()
            return output.split('\n')*.split(':').collect { it[0] }.join(", ")
        } catch (Exception e) {
            throw new Exception("Forbidden char check failed: $e - ${e.getCause()}")
        }
    }
}
