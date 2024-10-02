import FastaSequence

class FastaFile {
    static List<FastaSequence> parse(String fastaFilePath) {
        def sequences = []
        FastaSequence currentSequence = null
        new File(fastaFilePath).eachLine { line ->
            if (line.startsWith(">")) {
                if (currentSequence) {
                    sequences.add(currentSequence)
                }

                def header = line.substring(1).split(" ", 2).collect { it.trim() }
                currentSequence = new FastaSequence(
                    header[0], 
                    header.size() > 1 ? header[1] : "", 
                    ""
                )
            } else {
                currentSequence.sequence += line.trim()
            }
        }

        if (currentSequence) {
            sequences.add(currentSequence)
        }

        return sequences
    }

    static validate(String fastaFilePath, boolean isNucleic, Map appsConfig, List<String> appsToRun) {
        // Add application-specific forbidden characters
        String forbiddenChars = ""
        appsToRun.each{ app ->
            forbiddenChars += appsConfig[app].invalidChars ?: ""
        }

        // Remove duplicate characters
        forbiddenChars = forbiddenChars.toSet().join("")
        
        // Parse fasta file and check each sequence
        int numSequences = 0
        String errors = null
        FastaSequence currentSequence = null
        new File(fastaFilePath).eachLine { line ->
            if (line.startsWith(">")) {
                if (currentSequence) {
                    numSequences++
                    errors = currentSequence.validate(isNucleic, forbiddenChars)
                    if (errors) {
                        def errorMsg = "Forbidden characters found in sequence ${currentSequence.id}: ${errors}."
                        return [errorMsg, numSequences]
                    }
                }

                def header = line.substring(1).split(" ", 2).collect { it.trim() }
                currentSequence = new FastaSequence(
                    header[0], 
                    header.size() > 1 ? header[1] : "", 
                    ""
                )
            } else {
                currentSequence.sequence += line.trim()
            }
        }

        if (currentSequence) {
            numSequences++
            errors = currentSequence.validate(isNucleic, forbiddenChars)
            if (errors) {
                def errorMsg = "Forbidden characters found in sequence ${currentSequence.id}: ${errors}."
                return [errorMsg, numSequences]
            }
        }

        return [null, numSequences]
    }
}
