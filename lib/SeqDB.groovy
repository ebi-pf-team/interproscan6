import groovy.sql.Sql
import java.security.MessageDigest
import java.util.regex.Pattern

class SeqDB {
    private String path
    private Sql sql
    public static final Pattern eslDescription = ~/^source=(.+?)\s+coords=/

    SeqDB(String path) {
        this.path = path
        this.connect()
        this.createTables()
    }

   private void connect() {
        String url = "jdbc:sqlite:${this.path}"
        String driver = "org.sqlite.JDBC"
        try {
            this.sql = Sql.newInstance(url, driver)
        } catch (Exception e) {
            e.printStackTrace()
            throw e
        }
    }

    private void createTables() {
        String sql1 = """
            CREATE TABLE IF NOT EXISTS PROTEIN_SEQUENCE (
                md5          VARCHAR NOT NULL PRIMARY KEY,
                sequence     TEXT NOT NULL UNIQUE
            )
        """

        String sql2 = """
            CREATE TABLE IF NOT EXISTS PROTEIN (
                md5          VARCHAR NOT NULL,
                id           VARCHAR NOT NULL,
                description  VARCHAR NOT NULL,
                PRIMARY KEY (md5, id, description)
                FOREIGN KEY (md5) REFERENCES PROTEIN_SEQUENCE(md5)
            )
        """

        String sql3 = """
            CREATE TABLE IF NOT EXISTS NUCLEOTIDE_SEQUENCE (
                md5          VARCHAR NOT NULL PRIMARY KEY,
                sequence     TEXT NOT NULL
            )
        """

        String sql4 = """
            CREATE TABLE IF NOT EXISTS NUCLEOTIDE (
                md5          VARCHAR NOT NULL,
                id           VARCHAR NOT NULL,
                description  VARCHAR NOT NULL,
                PRIMARY KEY (md5, id, description),
                FOREIGN KEY (md5) REFERENCES NUCLEOTIDE_SEQUENCE(md5)
            )
        """

        String sql5 = "CREATE INDEX IF NOT EXISTS idx_nucleotide_id ON NUCLEOTIDE(id)"

        String sql6 = """
            CREATE TABLE IF NOT EXISTS PROTEIN_TO_NUCLEOTIDE (
                protein_md5  VARCHAR NOT NULL,
                nt_md5       VARCHAR NOT NULL,
                PRIMARY KEY (protein_md5, nt_md5)
                FOREIGN KEY (protein_md5) REFERENCES PROTEIN(md5),
                FOREIGN KEY (nt_md5) REFERENCES NUCLEOTIDE(md5) 
            )
        """

        this.sql.execute(sql1)
        this.sql.execute(sql2)
        this.sql.execute(sql3)
        this.sql.execute(sql4)
        this.sql.execute(sql5)
        this.sql.execute(sql6)
    }

    void close() {
        try {
            if (this.sql != null) {
                this.sql.close()
            }
        } catch (Exception e) {
            e.printStackTrace()
            throw e
        }
    }

    void loadFastaFile(String fastaFilePath, boolean isNucleic, boolean isTranslated) {
        File fastaFile = new File(fastaFilePath)
        def currentHeader = null
        def currentSeq = new StringBuilder()
        Map<String, Set<String>> ntSequences = [:]

        // Open a transaction
        sql.withTransaction {
            String table = isNucleic ? "NUCLEOTIDE" : "PROTEIN"
            String query1 = "INSERT OR IGNORE INTO ${table}_SEQUENCE (md5, sequence) VALUES (?, ?)"
            String query2 = "INSERT OR IGNORE INTO ${table} (md5, id, description) VALUES (?, ?, ?)"
            
            // Use batches for the sequence table
            sql.withBatch(100, query1) { ps1 ->
                // Use batches for the metadata table
                sql.withBatch(100, query2) { ps2 ->
                    // Read FASTA file and insert records
                    fastaFile.eachLine { line ->
                        line = line.trim()
                        if (line.startsWith(">")) {
                            if (currentHeader) {
                                String sequence = currentSeq.toString()
                                def seq = this.processRecord(currentHeader, sequence)
                                ps1.addBatch(seq.md5, sequence)
                                ps2.addBatch(seq.md5, seq.id, seq.description)

                                if (isTranslated) {
                                    // Extract the source nucleotide id and add to ntSequences
                                    def matcher = (seq.description =~ eslDescription)
                                    if (!matcher.find()) {
                                        throw new IllegalArgumentException("Invalid esl-translate FASTA header: ${seq.description}")
                                    }
                                    String sourceId = matcher.group(1)
                                    if (!ntSequences.containsKey(sourceId)) {
                                        ntSequences[sourceId] = [] as Set
                                    }
                                    ntSequences[sourceId] << seq.md5
                                }
                            }
                            currentHeader = line.substring(1).trim()
                            currentSeq = new StringBuilder()
                        } else {
                            currentSeq.append(line)
                        }
                    }
                    if (currentHeader) {
                        // Last record
                        String sequence = currentSeq.toString()
                        def seq = this.processRecord(currentHeader, sequence)
                        ps1.addBatch(seq.md5, sequence)
                        ps2.addBatch(seq.md5, seq.id, seq.description)

                        if (isTranslated) {
                            def matcher = (seq.description =~ eslDescription)
                            if (!matcher.find()) {
                                throw new IllegalArgumentException("Invalid esl-translate FASTA header: ${seq.description}")
                            }
                            String sourceId = matcher.group(1)
                            if (!ntSequences.containsKey(sourceId)) {
                                ntSequences[sourceId] = [] as Set
                            }
                            ntSequences[sourceId] << seq.md5
                        }
                    }
                }
            }

            if (isTranslated) {
                // Insert a mapping between the original nucleic sequences and their translation
                String query3 = "INSERT OR IGNORE INTO PROTEIN_TO_NUCLEOTIDE (protein_md5, nt_md5) VALUES (?, ?)"
                sql.withBatch(100, query3) { ps ->
                    // Get the MD5s of the original (nucleic) sequences using batches
                    ntSequences.keySet().toList().collate(100).each { List<String> chunk ->
                        String placeholders = chunk.collect { '?' }.join(',')
                        String query4 = "SELECT id, md5 FROM NUCLEOTIDE WHERE id IN (${placeholders})"
                        def records = []
                        sql.eachRow(query4, chunk) { row ->
                            String ntId = row.id
                            String ntMD5 = row.md5
                            ntSequences[ntId].each { protMD5 ->
                                records.add([protMD5, ntMD5])
                            } 
                        }

                        records.each {
                            ps.addBatch(it)
                        }
                    }
                }
            }
        }
    }

    private static Map processRecord(String header, String sequence) {
        def tokens = header.split(/\s+/, 2)
        String id = tokens[0]
        String description = tokens.length > 1 ? tokens[1] : ""

        MessageDigest md = MessageDigest.getInstance("MD5")
        byte[] digest = md.digest(sequence.toUpperCase().bytes)
        String md5 = digest.encodeHex().toString().toUpperCase()

        return [id: id, description: description, md5: md5]
    }

    void splitFasta(String outputPrefix, int maxSequencesPerFile) {
        String query = "SELECT md5, sequence FROM PROTEIN_SEQUENCE ORDER BY md5"
        int fileIndex = 1
        int seqCounter = 0

        // Open the first file for writing
        def writer = new File("${outputPrefix}.${fileIndex}.fasta").newWriter()

        this.sql.eachRow(query) { row ->
            if (seqCounter >= maxSequencesPerFile) {
                // We've reached the maximum for the current file: close it and open a new file
                writer.close()
                fileIndex++
                seqCounter = 0
                writer = new File("${outputPrefix}.${fileIndex}.fasta").newWriter()
            }

            writer.writeLine(">${row.md5}")

            String sequence = row.sequence
            for (int i = 0; i < sequence.length(); i += 60) {
                int end = Math.min(i + 60, sequence.length())
                writer.writeLine(sequence.substring(i, end))
            }
            seqCounter++
        }

        writer.close()
    }
}