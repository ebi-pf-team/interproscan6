import groovy.sql.Sql
import java.security.MessageDigest
import java.util.regex.Pattern

class SeqDB {
    // This class contains the methods for querying the internal IPS6 seq db, and handling any errors
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

    void splitFasta(String outputPrefix, int maxSequencesPerFile, boolean nucleic) {
        String query = nucleic ? """SELECT P2N.nt_md5, S.md5, S.sequence
            FROM PROTEIN_SEQUENCE AS S
            INNER JOIN PROTEIN AS P ON S.md5 = P.md5
            INNER JOIN PROTEIN_TO_NUCLEOTIDE AS P2N ON P.md5 = P2N.protein_md5
            ORDER BY P2N.nt_md5""" : "SELECT NULL AS nt_md5, md5, sequence FROM PROTEIN_SEQUENCE ORDER BY md5"
        int fileIndex = 1
        def batch = []
        def currentMD5 = null
        def seenProteinMd5s = [] as Set
        def writer = null

        this.sql.eachRow(query) { row ->
            // If we encounter a new nt_md5 and the batch is full --> write to a new batch
            if (currentMD5 && currentMD5 != row.nt_md5 && batch.size() >= maxSequencesPerFile) {
                writer = new File("${outputPrefix}.${fileIndex}.fasta").newWriter()
                for (record in batch) {
                    writer.writeLine(">${record.md5}")
                    String sequence = record.sequence
                    for (int i = 0; i < sequence.length(); i += 60) {
                        int end = Math.min(i + 60, sequence.length())
                        writer.writeLine(sequence.substring(i, end))
                    }
                }
                writer.close()
                fileIndex++
                batch.clear()
                seenProteinMd5s.clear()
            }

            if (!seenProteinMd5s.contains(row.md5)) {
                seenProteinMd5s.add(row.md5)
                // Don't add the GroovyResultSet item, else when writing the final batch we will get a ResultSet closed error
                batch.add([
                    md5     : row.md5,
                    sequence: row.sequence
                ])
            }

            currentMD5 = nucleic ? row.nt_md5 : row.md5
        }

        // write the final batch
        if (!batch.isEmpty()) {
            writer = new File("${outputPrefix}.${fileIndex}.fasta").newWriter()
            batch.each { record ->
                writer.writeLine(">${record.md5}")
                def seq = record.sequence
                for (int i = 0; i < seq.length(); i += 60) {
                    writer.writeLine(seq.substring(i, Math.min(i + 60, seq.length())))
                }
            }
            writer.close()
        }
    }

    def groupProteins(Map proteinMatches) {
        /* Gather nucleotide Seq IDs and child protein MD5s so we can gather all ORFs from the same
        parent NT seq together in the final output. */
        def nucleicRelationships = [:]  // [ntSeqId: [proteinMd5]]
        proteinMatches.each { proteinMD5, matchesNode ->
            // get all parent NT seq Ids
            def query = """SELECT DISTINCT nt_md5
                FROM PROTEIN_TO_NUCLEOTIDE
                WHERE protein_md5 = ?
                """
            this.sql.eachRow(query, [proteinMD5]) { row ->
                nucleicRelationships.computeIfAbsent(row.nt_md5, { [] as Set })
                nucleicRelationships[row.nt_md5].add(proteinMD5)
            }
        }
        return nucleicRelationships
    }

    List<String> proteinMd5ToNucleicSeq(String proteinMD5) {
        def query = """SELECT N.id AS nid, N.description, S.sequence, P.id AS pid
            FROM NUCLEOTIDE AS N
            INNER JOIN NUCLEOTIDE_SEQUENCE AS S ON N.md5 = S.md5
            INNER JOIN PROTEIN_TO_NUCLEOTIDE AS N2P ON N.md5 = N2P.nt_md5
            INNER JOIN PROTEIN AS P ON N2P.protein_md5 = P.md5
            WHERE P.md5 = ?;
            """
        return this.sql.rows(query, [proteinMD5])
    }

    List<String> proteinMd5ToProteinSeq(String proteinMD5) {
        def query = """SELECT P.id, P.description, S.sequence
            FROM PROTEIN AS P
            INNER JOIN PROTEIN_SEQUENCE AS S ON P.md5 = S.md5
            WHERE P.md5 = ?;
            """
        return this.sql.rows(query, [proteinMD5])
    }

    List<String> nucleicMd5ToNucleicSeq(String nucleicMD5) {
        def query = """SELECT N.id, N.description, S.sequence
            FROM NUCLEOTIDE AS N
            INNER JOIN NUCLEOTIDE_SEQUENCE AS S ON N.md5 = S.md5
            WHERE N.md5 = ?;
            """
        return this.sql.rows(query, [nucleicMD5])
    }

    List<String> getOrfSeq(String proteinMD5, String nucleicMD5) {
        /* A protein seq may be associated with multiple nucleotide seqs.
        This could because of duplication or multiple nucleotide seqs encoding the same protein.
        To ensure we return the data for the protein associated for ONLY the current working nucleotide
        seq we need to check the nucleotide seq ID against the description generated by ESL_translate.
         */
        def query = """SELECT P.id, P.description, S.sequence, N.id AS nt_id
            FROM PROTEIN AS P
            INNER JOIN PROTEIN_SEQUENCE AS S ON P.md5 = S.md5
            INNER JOIN PROTEIN_TO_NUCLEOTIDE AS N2P ON P.md5 = N2P.protein_md5
            INNER JOIN NUCLEOTIDE AS N ON N2P.nt_md5 = N.md5
            WHERE P.md5 = ? AND N.md5 = ? AND P.description LIKE 'source=' || N.id || '%';
            """
        return this.sql.rows(query, [proteinMD5, nucleicMD5])
    }
}