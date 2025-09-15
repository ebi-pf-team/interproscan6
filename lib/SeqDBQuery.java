import java.sql.*;
import java.util.*;
import java.security.MessageDigest;
java.util.regex.Pattern;

public class SeqDB {
    private Connection connection;
    private static final int insertBatchSize = 100;
    private static final Pattern eslDescriptionPattern = Pattern.compile("^source=(.+?)\\s+coords=");

    public SeqDB(String dbPath) throws SQLException {
        try {
            Class.forName("org.sqlite.JDBC");
        } catch (ClassNotFoundException e) {
            throw new SQLException("SQLite JDBC driver not found", e);
        }
        
        String url = "jdbc:sqlite:" + dbPath;
        this.connection = DriverManager.getConnection(url);
        
        // Enable (potential) performance optimizations
        connection.createStatement().execute("PRAGMA cache_size = 10000");
        connection.createStatement().execute("PRAGMA temp_store = MEMORY");
        connection.createStatement().execute("PRAGMA synchronous = OFF");
        connection.createStatement().execute("PRAGMA journal_mode = MEMORY");
    }

    public void createTables() throws SQLException {
        String createProteinTable = """
            CREATE TABLE IF NOT EXISTS PROTEIN (
                id TEXT PRIMARY KEY,
                md5 TEXT,
                description TEXT
            );
            """;

        String createProteinSeqTable = """
            CREATE TABLE IF NOT EXISTS PROTEIN_SEQUENCE (
                md5 TEXT PRIMARY KEY,
                sequence TEXT
            );
            """;

        String createNucleotideTable = """
            CREATE TABLE IF NOT EXISTS NUCLEOTIDE (
                id TEXT PRIMARY KEY,
                md5 TEXT
            );
            """;

        String createNucleotideSeqTable = """
            CREATE TABLE IF NOT EXISTS NUCLEOTIDE_SEQUENCE (
                md5 TEXT PRIMARY KEY,
                sequence TEXT
            );
            """;

        String createProteinToNucleotideTable = """
            CREATE TABLE IF NOT EXISTS PROTEIN_TO_NUCLEOTIDE (
                protein_md5 TEXT,
                nt_md5 TEXT,
                PRIMARY KEY (protein_md5, nt_md5)
            );
            """;

        try (Statement stmt = connection.createStatement()) {
            stmt.execute(createProteinTable);
            stmt.execute(createProteinSeqTable);
            stmt.execute(createNucleotideTable);
            stmt.execute(createNucleotideSeqTable);
            stmt.execute(createProteinToNucleotideTable);
        }
    }


    public void loadFastaFile(String fastaFilePath, boolean isNucleic, boolean isTranslated) {
        File fastaFile = new File(fastaFilePath);
        String currentHeaeder = null;
        String currentSeq = null;
        Map<String, Set<String>> ntSequences = [:];

        String table = isNucleic ? "NUCLEOTIDE" : "PROTEIN";
        String insertSeqQuery = "INSERT OR IGNORE INTO " + table + " (id, md5) VALUES (?, ?)";
        String insertSeqDataQuery = "INSERT OR IGNORE INTO " + table + "(md5, id, description) VALUES (?, ?, ?)";

        int recordCount = 0;

        try (PreparedStatement stmt1 = connection.prepareStatement(insertSeqQuery)) {
            try (PreparedStatement stmt2 = connection.prepareStatement(insertSeqDataQuery)) {
                try (BufferedReader br = new BufferedReader(new FileReader(fastaFile))) {
                    String line;
                    while ((line = br.readLine()) != null) {
                        if (line.startsWith(">")) {
                            if (currentHeader) { // Process the previous record
                                String sequence = currentSeq.toString();
                                Map<String, String> seqData = this.processRecord(currentHeader, sequence)
                                stmt1.setString(1, seqData.get("md5"));
                                stmt1.setString(2, seqData.get("sequence"));
                                stmt1.addBatch();
                                stmt2.setString(1, seqData.get("md5"));
                                stmt2.setString(2, seqData.get("id"));
                                stmt2.setString(3, seqData.get("description"));
                                stmt2.addBatch();

                                if (isTranslated)  {
                                    // Extract the source nucleotide id and add to ntSequences
                                    Matcher matcher = eslDescriptionPattern.matcher(seqData.get("description"));
                                    if (matcher.find()) {
                                        String sourceId = matcher.group(1);
                                        ntSequences.computeIfAbsent(sourceId, k -> new HashSet<>()).add(seqData.get("md5"));
                                    } else {
                                        throw new IllegalArgumentException("Invalid esl-translate FASTA header: " + seqData.get("description"));
                                    }
                                }

                                recordCount++;
                                if (recordCount % insertBatchSize == 0) {
                                    stmt1.executeBatch();
                                    stmt2.executeBatch();
                                    stmt1.clearBatch();
                                    stmt2.clearBatch();
                                    recordCount = 0;
                                }
                            }
                            currentHeader = line.substring(1).trim();
                            currentSeq = new StringBuilder();
                        }
                    }
                } // end of buffered reader

                // Finished parsing the file, now process the last recrd
                if (currentHeader) {
                    String sequence = currentSeq.toString();
                    Map<String, String> seqData = this.processRecord(currentHeader, sequence)
                    stmt1.setString(1, seqData.get("md5"));
                    stmt1.setString(2, seqData.get("sequence"));
                    stmt1.addBatch();
                    stmt2.setString(1, seqData.get("md5"));
                    stmt2.setString(2, seqData.get("id"));
                    stmt2.setString(3, seqData.get("description"));
                    stmt2.addBatch();

                    if (isTranslated)  {
                        // Extract the source nucleotide id and add to ntSequences
                        Matcher matcher = eslDescriptionPattern.matcher(seqData.get("description"));
                        if (matcher.find()) {
                            String sourceId = matcher.group(1);
                            ntSequences.computeIfAbsent(sourceId, k -> new HashSet<>()).add(seqData.get("md5"));
                        } else {
                            throw new IllegalArgumentException("Invalid esl-translate FASTA header: " + seqData.get("description"));
                        }
                    }
                }

                // insert any remaining batched records
                stmt1.executeBatch();
                stmt2.executeBatch();
                stmt1.clearBatch();
                stmt2.clearBatch();

            } // end of prepared statement 2
        } // end of prepared statement 1


    }

    private static Map processRecord(String header, String sequences) {
        List<String> tokens = header.split("\\s+", 2);
        String id = tokens[0];
        String description = tokens.length > 1 ? tokens[1] : "";
        MessageDigest md = MessageDigest.getInstance("MD5");
        byte[] digest = md.digest(sequence.toUpperCase().getBytes(StandardCharsets.UTF_8));
        String md5 = bytesToHex(digest);
        return Map.of("id", id, "description", description, "md5", md5, "sequence", sequence);
    }

    public Map<String, List<ProteinData>> proteinMd5sToProteinSeqs(List<String> proteinMD5s) throws SQLException {
        StringBuilder placeholders = new StringBuilder();
        for (int i = 0; i < proteinMD5s.size(); i++) {
            if (i > 0) placeholders.append(',');
            placeholders.append('?');
        }

        String query = String.format(
            """
            SELECT P.id, P.description, S.sequence, P.md5
            FROM PROTEIN AS P
            INNER JOIN PROTEIN_SEQUENCE AS S ON P.md5 = S.md5
            WHERE P.md5 IN (%s)
            """, 
            placeholders
        );

        Map<String, List<ProteinData>> resultMap = new HashMap<>(proteinMD5s.size());

        try (PreparedStatement ps = connection.prepareStatement(query)) {
            for (int i = 0; i < proteinMD5s.size(); i++) {
                ps.setString(i + 1, proteinMD5s.get(i));
            }

            try (ResultSet rs = ps.executeQuery()) {
                while (rs.next()) {
                    String md5 = rs.getString("md5");
                    ProteinData protein = new ProteinData(
                        rs.getString("id"),
                        rs.getString("description"), 
                        rs.getString("sequence")
                    );

                    resultMap.computeIfAbsent(md5, k -> new ArrayList<>()).add(protein);
                }
            }
        }

        return resultMap;
    }

    public List<NucleotideData> nucleicMd5ToNucleicSeq(String nucleicMd5) throws SQLException {
        List<NucleotideData> result = new ArrayList<>();
        
        String sql = """
            SELECT N.id, S.sequence
            FROM NUCLEOTIDE AS N
            INNER JOIN NUCLEOTIDE_SEQUENCE AS S ON N.md5 = S.md5
            WHERE N.md5 = ?
            """;
        
        try (PreparedStatement stmt = connection.prepareStatement(sql)) {
            stmt.setString(1, nucleicMd5);

            try (ResultSet rs = stmt.executeQuery()) {
                while (rs.next()) {
                    String id = rs.getString("id");
                    String sequence = rs.getString("sequence");
                    result.add(new NucleotideData(id, sequence));
                }
            }
        }

        return result;
    }

    public Map<String, List<NucleotideData>> nucleicMd5sToNucleicSeqs(List<String> nucleicMd5s) throws SQLException {
        if (nucleicMd5s == null || nucleicMd5s.isEmpty()) {
            return new HashMap<>();
        }

        Map<String, List<NucleotideData>> result = new HashMap<>();
        
        StringBuilder placeholders = new StringBuilder();
        for (int i = 0; i < nucleicMd5s.size(); i++) {
            if (i > 0) placeholders.append(",");
            placeholders.append("?");
        }

        String sql = String.format("""
            SELECT N.md5, N.id, S.sequence
            FROM NUCLEOTIDE AS N
            INNER JOIN NUCLEOTIDE_SEQUENCE AS S ON N.md5 = S.md5
            WHERE N.md5 IN (%s)
            """, placeholders);
        
        try (PreparedStatement stmt = connection.prepareStatement(sql)) {
            for (int i = 0; i < nucleicMd5s.size(); i++) {
                stmt.setString(i + 1, nucleicMd5s.get(i));
            }

            try (ResultSet rs = stmt.executeQuery()) {
                while (rs.next()) {
                    String md5 = rs.getString("md5");
                    String id = rs.getString("id");
                    String sequence = rs.getString("sequence");

                    result.computeIfAbsent(md5, k -> new ArrayList<>())
                          .add(new NucleotideData(id, sequence));
                }
            }
        }

        return result;
    }

    public List<ProteinData> getOrfSeq(String proteinMd5, String nucleicMd5) throws SQLException {
        List<ProteinData> result = new ArrayList<>();
        
        String sql = """
            SELECT P.id, P.description, S.sequence, N.id AS nt_id
            FROM PROTEIN AS P
            INNER JOIN PROTEIN_SEQUENCE AS S ON P.md5 = S.md5
            INNER JOIN PROTEIN_TO_NUCLEOTIDE AS N2P ON P.md5 = N2P.protein_md5
            INNER JOIN NUCLEOTIDE AS N ON N2P.nt_md5 = N.md5
            WHERE P.md5 = ? AND N.md5 = ? AND P.description LIKE 'source=' || N.id || '%'
            """;
        
        try (PreparedStatement stmt = connection.prepareStatement(sql)) {
            stmt.setString(1, proteinMd5);
            stmt.setString(2, nucleicMd5);

            try (ResultSet rs = stmt.executeQuery()) {
                while (rs.next()) {
                    String id = rs.getString("id");
                    String description = rs.getString("description");
                    String sequence = rs.getString("sequence");
                    result.add(new ProteinData(id, description, sequence));
                }
            }
        }

        return result;
    }

    public Map<String, Map<String, List<ProteinData>>> getOrfSeqsBulk(List<String> proteinMd5s, List<String> nucleicMd5s) throws SQLException {
        if (proteinMd5s == null || proteinMd5s.isEmpty() || nucleicMd5s == null || nucleicMd5s.isEmpty()) {
            return new HashMap<>();
        }

        Map<String, Map<String, List<ProteinData>>> result = new HashMap<>();
        
        StringBuilder proteinPlaceholders = new StringBuilder();
        for (int i = 0; i < proteinMd5s.size(); i++) {
            if (i > 0) proteinPlaceholders.append(",");
            proteinPlaceholders.append("?");
        }

        StringBuilder nucleicPlaceholders = new StringBuilder();
        for (int i = 0; i < nucleicMd5s.size(); i++) {
            if (i > 0) nucleicPlaceholders.append(",");
            nucleicPlaceholders.append("?");
        }

        String sql = String.format("""
            SELECT P.md5 as protein_md5, N.md5 as nucleic_md5, P.id, P.description, S.sequence
            FROM PROTEIN AS P
            INNER JOIN PROTEIN_SEQUENCE AS S ON P.md5 = S.md5
            INNER JOIN PROTEIN_TO_NUCLEOTIDE AS N2P ON P.md5 = N2P.protein_md5
            INNER JOIN NUCLEOTIDE AS N ON N2P.nt_md5 = N.md5
            WHERE P.md5 IN (%s) AND N.md5 IN (%s) AND P.description LIKE 'source=' || N.id || '%%'
            """, proteinPlaceholders, nucleicPlaceholders);
        
        try (PreparedStatement stmt = connection.prepareStatement(sql)) {
            int paramIndex = 1;
            for (String proteinMd5 : proteinMd5s) {
                stmt.setString(paramIndex++, proteinMd5);
            }
            
            for (String nucleicMd5 : nucleicMd5s) {
                stmt.setString(paramIndex++, nucleicMd5);
            }

            try (ResultSet rs = stmt.executeQuery()) {
                while (rs.next()) {
                    String proteinMd5 = rs.getString("protein_md5");
                    String nucleicMd5 = rs.getString("nucleic_md5");
                    String id = rs.getString("id");
                    String description = rs.getString("description");
                    String sequence = rs.getString("sequence");

                    result.computeIfAbsent(proteinMd5, k -> new HashMap<>())
                          .computeIfAbsent(nucleicMd5, k -> new ArrayList<>())
                          .add(new ProteinData(id, description, sequence));
                }
            }
        }

        return result;
    }
    public Map<String, Set<String>> groupProteins(Map<String, Object> proteins) throws SQLException {
        Map<String, Set<String>> nucleicToProteins = new HashMap<>();
        Set<String> proteinMd5s = proteins.keySet();
        
        if (proteinMd5s.isEmpty()) {
            return nucleicToProteins;
        }
        
        StringBuilder placeholders = new StringBuilder();
        List<String> proteinMd5List = new ArrayList<>(proteinMd5s);
        for (int i = 0; i < proteinMd5List.size(); i++) {
            if (i > 0) placeholders.append(",");
            placeholders.append("?");
        }

        String sql = String.format("""
            SELECT DISTINCT protein_md5, nt_md5
            FROM PROTEIN_TO_NUCLEOTIDE
            WHERE protein_md5 IN (%s)
            """, placeholders);
        
        try (PreparedStatement stmt = connection.prepareStatement(sql)) {
            for (int i = 0; i < proteinMd5List.size(); i++) {
                stmt.setString(i + 1, proteinMd5List.get(i));
            }

            try (ResultSet rs = stmt.executeQuery()) {
                while (rs.next()) {
                    String proteinMd5 = rs.getString("protein_md5");
                    String nucleicMd5 = rs.getString("nt_md5");

                    nucleicToProteins.computeIfAbsent(nucleicMd5, k -> new HashSet<>())
                                   .add(proteinMd5);
                }
            }
        }

        return nucleicToProteins;
    }

    public Map<String, Set<String>> groupProteinsBulk(Set<String> allProteinMd5s) throws SQLException {
        if (allProteinMd5s == null || allProteinMd5s.isEmpty()) {
            return new HashMap<>();
        }

        Map<String, Set<String>> nucleicToProteins = new HashMap<>();
        List<String> proteinMd5List = new ArrayList<>(allProteinMd5s);
        
        StringBuilder placeholders = new StringBuilder();
        for (int i = 0; i < proteinMd5List.size(); i++) {
            if (i > 0) placeholders.append(",");
            placeholders.append("?");
        }

        String sql = String.format("""
            SELECT DISTINCT protein_md5, nt_md5
            FROM PROTEIN_TO_NUCLEOTIDE
            WHERE protein_md5 IN (%s)
            """, placeholders);
        
        try (PreparedStatement stmt = connection.prepareStatement(sql)) {
            for (int i = 0; i < proteinMd5List.size(); i++) {
                stmt.setString(i + 1, proteinMd5List.get(i));
            }

            try (ResultSet rs = stmt.executeQuery()) {
                while (rs.next()) {
                    String proteinMd5 = rs.getString("protein_md5");
                    String nucleicMd5 = rs.getString("nt_md5");

                    nucleicToProteins.computeIfAbsent(nucleicMd5, k -> new HashSet<>())
                                   .add(proteinMd5);
                }
            }
        }

        return nucleicToProteins;
    }

    public void close() throws SQLException {
        if (connection != null && !connection.isClosed()) {
            connection.close();
        }
    }
}

class NucleotideData {
    public final String id;
    public final String sequence;

    public NucleotideData(String id, String sequence) {
        this.id = id;
        this.sequence = sequence;
    }
}

record ProteinData(String id, String description, String sequence) {}
