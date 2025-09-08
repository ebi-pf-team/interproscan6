import java.sql.*;
import java.util.*;

public class SeqDBQuery {
    private Connection connection;

    public SeqDBQuery(String dbPath) throws SQLException {
        try {
            Class.forName("org.sqlite.JDBC");
        } catch (ClassNotFoundException e) {
            throw new SQLException("SQLite JDBC driver not found", e);
        }
        
        String url = "jdbc:sqlite:" + dbPath;
        this.connection = DriverManager.getConnection(url);
        
        // Enable performance optimizations
        connection.createStatement().execute("PRAGMA cache_size = 10000");
        connection.createStatement().execute("PRAGMA temp_store = MEMORY");
        connection.createStatement().execute("PRAGMA synchronous = OFF");
        connection.createStatement().execute("PRAGMA journal_mode = MEMORY");
    }

    // Fixed to match Groovy SeqDB.proteinMd5sToProteinSeqs()
    public Map<String, List<ProteinData>> proteinMd5sToProteinSeqs(List<String> proteinMD5s) throws SQLException {
        if (proteinMD5s == null || proteinMD5s.isEmpty()) {
            return new HashMap<>();
        }

        StringBuilder placeholders = new StringBuilder();
        for (int i = 0; i < proteinMD5s.size(); i++) {
            if (i > 0) placeholders.append(',');
            placeholders.append('?');
        }

        // Fixed table names to match your Groovy schema
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

    // Fixed to match Groovy SeqDB.nucleicMd5ToNucleicSeq()
    public List<NucleotideData> nucleicMd5ToNucleicSeq(String nucleicMd5) throws SQLException {
        List<NucleotideData> result = new ArrayList<>();
        
        // Fixed table names and query structure
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

    // Bulk version of nucleicMd5ToNucleicSeq for performance
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

    // Fixed to match Groovy SeqDB.getOrfSeq()
    public List<ProteinData> getOrfSeq(String proteinMd5, String nucleicMd5) throws SQLException {
        List<ProteinData> result = new ArrayList<>();
        
        // Fixed query to match Groovy version exactly
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

    // Bulk version of getOrfSeq for performance
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

        // Fixed query to match Groovy getOrfSeq logic
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

    // Fixed to match Groovy SeqDB.groupProteins()
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

        // Fixed query to match Groovy groupProteins exactly
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

    // Bulk version of groupProteins for performance
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
