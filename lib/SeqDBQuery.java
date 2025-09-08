import java.sql.*;
import java.util.*;

public class SeqDBQuery {
    private Connection connection;
    // Remove this line - it's not used and has a typo
    // private PrepareStatement bulkProteinQuery;

    public SeqDBQuery(String dbPath) throws SQLException {
        try {
            // Explicitly load the SQLite JDBC driver
            Class.forName("org.sqlite.JDBC");
        } catch (ClassNotFoundException e) {
            throw new SQLException("SQLite JDBC driver not found", e);
        }

        String url = "jdbc:sqlite:" + dbPath;
        this.connection = DriverManager.getConnection(url);
        // Enable performance optimizations
        connection.createStatement().execute("PRAGMA cache_size = 10000");
        connection.createStatement().execute("PRAGMA temp_store = MEMORY");
    }

    // Fix: change "string" to "String" (capital S)
    public Map<String, List<ProteinData>> proteinMd5sToProteinSeqs(List<String> proteinMD5s) throws SQLException {
        if (proteinMD5s == null || proteinMD5s.isEmpty()) {
            return new HashMap<>();
        }

        StringBuilder placeholders = new StringBuilder();
        for (int i = 0; i < proteinMD5s.size(); i++) {  // Fix: add space before ")"
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
            // Set parameters
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
    
    public void close() throws SQLException {
        if (connection != null && !connection.isClosed()) {
            connection.close();
        }
    }
}

// Define the record outside the class or at the end
record ProteinData(String id, String description, String sequence) {}
