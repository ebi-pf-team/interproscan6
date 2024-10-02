import java.security.MessageDigest

class FastaSequence {
    // DNA, RNA, gaps
    static final String NUCLEIC_ALPHABET = ~/(?i)[ACGTUN\-\.\s]+/
    // 20 standard AAs, Sec, Pyl, any/unknown, Asx, Glx, Xle, gaps
    static final String PROTEIN_ALPHABET = /(?i)[ACDEFGHIKLMNPQRSTVWYUOXBZJ\-\.\s]+/

    String id
    String description
    String sequence

    String getMD5() {
        MessageDigest md = MessageDigest.getInstance("MD5")
        byte[] digest = md.digest(sequence.bytes)
        return digest.encodeHex().toString()
    }

    String validate(boolean isNucleic = false, String extraForbiddenChars = "") {
        def alphabet = isNucleic ? this.NUCLEIC_ALPHABET : this.PROTEIN_ALPHABET
        def pattern = "[^${alphabet}]"
        def invalidChars = this.sequence.findAll(/(?i)${pattern}/).unique().join("")
        if (invalidChars) {
            return invalidChars
        }

        if (!isNucleic && extraForbiddenChars) {
            pattern = /(?i)[\Q${extraForbiddenChars}\E]/
            invalidChars = this.sequence.findAll(pattern).unique().join("")
            if (invalidChars) {
                return invalidChars
            }
        }
        
        return null
    }

    static void checkFastaFile(String fastaFilePath, boolean isNucleic, Map appsConfig, List<String> appsToRun) {
        // Add application-specific forbidden characters
        String forbiddenChars = ""
        appsToRun.each{ app ->
            forbiddenChars += appsConfig[app].invalidChars ?: ""
        }

        // Remove duplicate characters
        forbiddenChars = forbiddenChars.toSet().join("")
        
        // Parse fasta file and check each sequence
        String errors = null
        FastaSequence currentSequence = null
        new File(fastaFilePath).eachLine { line ->
            if (line.startsWith(">")) {
                if (currentSequence) {
                    errors = currentSequence.validate(isNucleic, forbiddenChars)
                    if (errors) {
                        println "Error: forbidden characters found in sequence ${currentSequence.id}: ${errors}"
                        System.exit(1)
                    }
                }

                def header = line.substring(1).split(" ", 2).collect { it.trim() }
                currentSequence = new FastaSequence(
                    id: header[0], 
                    description: header.size() > 1 ? header[1] : "", 
                    sequence: ""
                )
            } else {
                currentSequence.sequence += line.trim()
            }
        }

        if (currentSequence) {
            errors = currentSequence.validate(isNucleic, forbiddenChars)
            if (errors) {
                println "Error: forbidden characters found in sequence ${currentSequence.id}: ${errors}"
                System.exit(1)
            }
        }
    }
}
