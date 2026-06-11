package uk.ac.ebi.interpro
import java.nio.file.Path

class FastaFile {
    // DNA, RNA, gaps
    static final String NUCLEIC_ALPHABET = "ACGTURYMKSWHBVDN\\.\\s"
    // 20 standard AAs, Sec, Pyl, any/unknown, Asx, Glx, Xle
    static final String PROTEIN_ALPHABET = "ACDEFGHIKLMNPQRSTVWYUOXBZJ\\s"

    static Map<String, String> parse(Path fastaFile) {
        def sequences = [:]
        def md5 = null
        StringBuilder builder = null

        fastaFile.eachLine { line ->
            if (line.startsWith(">")) {
                if (builder) {
                    sequences[md5] = builder.toString()
                }
                md5 = line.substring(1).trim()
                builder = new StringBuilder()
            } else {
                builder.append(line.replaceAll("\\s+", ""))
            }
        }
        if (builder) {
            sequences[md5] = builder.toString()
        }
        return sequences
    }

    static String validate(Path fastaFile, boolean isNucleic) {
        String allowed = isNucleic ? NUCLEIC_ALPHABET : PROTEIN_ALPHABET
        def pattern = ~"^[${allowed}]+\$"
        BufferedReader reader = fastaFile.newReader()
        String seqId
        String line

        try {
            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {
                    seqId = line.substring(1).trim()
                } else {
                    String seq = line.replaceAll("\\s+", "")
                    if (seq.isEmpty()) {
                        continue
                    } else if (!pattern.matcher(seq).matches()) {
                        return seqId
                    }
                }
            }
        } finally {
            reader.close()
        }

        return null
    }

    static void write(Path inputFasta,
                      Path outputFasta,
                      String extraForbidden,
                      Map<Character,Character> substitutions) {
        Set<Character> baseAllowed = PROTEIN_ALPHABET.toList() as Set
        Set<Character> forbid = extraForbidden.toList() as Set

        inputFasta.withReader { reader ->
            outputFasta.withWriter { writer ->
                String header = null
                StringBuilder builder = new StringBuilder()

                reader.eachLine { line ->
                    if (line.startsWith(">")) {
                        if (header) {
                            processEntry(header, builder, baseAllowed, forbid, substitutions, writer)
                        }
                        header = line.trim()
                        builder = new StringBuilder()
                    } else {
                        builder.append(line.replaceAll("\\s+", ""))
                    }
                }
                if (header) {
                    processEntry(header, builder, baseAllowed, forbid, substitutions, writer)
                }
            }
        }
    }

    private static void processEntry(String header,
                                     StringBuilder inSeq,
                                     Set<Character> baseAllowed,
                                     Set<Character> forbid,
                                     Map<Character,Character> substitutions,
                                     Writer writer) {
        StringBuilder outSeq = new StringBuilder()
        boolean ok = true

        inSeq.toString().each { ch ->
            // substitute if needed
            if (substitutions.containsKey(ch)) {
                ch = substitutions[ch]
            }
            // check forbidden or outside base alphabet
            if (forbid.contains(ch) || !baseAllowed.contains(ch)) {
                ok = false
                return
            }
            outSeq.append(ch)
        }

        if (ok) {
            writer.writeLine(header)
            outSeq.toString()
                  .eachMatch(/.{1,60}/) { writer.writeLine(it) }
        }
    }

    static List<Map> chunkSequences(Path path, int chunkSize, int overlap) {
        // Parse sequences
        def sequences = this.parse(path)

        // Chunk sequences
        def result = sequences.collect { seqId, sequence ->
            def seqLength = sequence.length()

            if (chunkSize <= 0 || seqLength <= chunkSize) {
                return [id: seqId, chunks: [sequence]]
            }

            def chunks = []
            int step = chunkSize - overlap
            for (int start = 0; start < seqLength; start += step) {
                int end = Math.min(start + chunkSize, seqLength)
                chunks << sequence.substring(start, end)
                if (end == seqLength) break
            }

            return [id: seqId, chunks: chunks]
        }

        return result
    }
}
