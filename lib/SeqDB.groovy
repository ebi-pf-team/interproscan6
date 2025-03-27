import groovy.sql.Sql
import java.security.MessageDigest

class SeqDB {
    private String path
    private Sql sql

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
                protein_md5  VARCHAR PRIMARY KEY,
                sequence     TEXT UNIQUE
            )
        """

        String sql2 = """
            CREATE TABLE IF NOT EXISTS PROTEIN (
                protein_md5  VARCHAR,
                id           VARCHAR,
                description  VARCHAR,
                FOREIGN KEY (protein_md5) REFERENCES PROTEIN_SEQUENCE(protein_md5),
                UNIQUE (protein_md5, id, description)
            )
        """

        String sql3 = """
            CREATE TABLE IF NOT EXISTS NUCLEOTIDE_SEQUENCE (
                nt_md5       VARCHAR PRIMARY KEY,
                sequence     TEXT
            )
        """

        String sql4 = """
            CREATE TABLE IF NOT EXISTS NUCLEOTIDE (
                nt_md5       VARCHAR,
                id           VARCHAR,
                description  VARCHAR,
                FOREIGN KEY (nt_md5) REFERENCES NUCLEOTIDE_SEQUENCE(nt_md5),
                UNIQUE (nt_md5, id, description)
            )
        """

        String sql5 = """
            CREATE TABLE IF NOT EXISTS PROTEIN_TO_NUCLEOTIDE (
                protein_md5 VARCHAR,
                nt_md5 VARCHAR,
                FOREIGN KEY (protein_md5) REFERENCES PROTEIN(protein_md5),
                FOREIGN KEY (nt_md5) REFERENCES NUCLEOTIDE(nt_md5),
                UNIQUE (protein_md5, nt_md5)
            )
        """

        this.sql.execute(sql1)
        this.sql.execute(sql2)
        this.sql.execute(sql3)
        this.sql.execute(sql4)
        this.sql.execute(sql5)
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

    void loadFastaFile(String fastaFilePath, boolean isNucleotide) {
        File fastaFile = new File(fastaFilePath)
        def currentHeader = null
        def currentSeq = new StringBuilder()
        fastaFile.eachLine { line ->
            line = line.trim()
            if (line.startsWith(">")) {
                if (currentHeader) {
                    this.insertSequence(currentHeader, currentSeq.toString(), isNucleotide)
                }
                currentHeader = line.substring(1).trim()
                currentSeq = new StringBuilder()
            } else {
                currentSeq.append(line)
            }
        }
        if (currentHeader) {
            this.insertSequence(currentHeader, currentSeq.toString(), isNucleotide)
        }
    }

    private insertSequence(String header, String sequence, boolean isNucleotide) {
        // Assume the header is formatted as "id description..."
        def tokens = header.split(/\s+/, 2)
        String id = tokens[0]
        String description = tokens.length > 1 ? tokens[1] : ""
        String md5 = this.computeMD5(sequence)

        if (isNucleotide) {
            this.sql.execute("INSERT OR IGNORE INTO NUCLEOTIDE_SEQUENCE (nt_md5, sequence) VALUES (?, ?)", [md5, sequence])
            this.sql.execute("INSERT OR IGNORE INTO NUCLEOTIDE (nt_md5, id, description) VALUES (?, ?, ?)", [md5, id, description])
        } else {
            this.sql.execute("INSERT OR IGNORE INTO PROTEIN_SEQUENCE (protein_md5, sequence) VALUES (?, ?)", [md5, sequence])
            this.sql.execute("INSERT OR IGNORE INTO PROTEIN (protein_md5, id, description) VALUES (?, ?, ?)", [md5, id, description])
        }
    }

    private static String computeMD5(String sequence) {
        MessageDigest md = MessageDigest.getInstance("MD5")
        byte[] digest = md.digest(sequence.toUpperCase().bytes)
        return digest.encodeHex().toString().toUpperCase()
    }
}