class SeqDatabase {
    // This class contains the method for querying the internal IPS6 seq db, and handling any errors
    String dbPath

    SeqDatabase(String dbPath) {
        this.dbPath = dbPath
    }

    static List<String> query(String query) {
        try {
            def cmd = ["sqlite3", "--tabs", this.dbPath, query]
            def process = cmd.execute()
            process.waitFor()
            if (process.exitValue() != 0) {
                throw new Exception("SQLite query failed: ${process.err.text.trim()}")
            }
            def output = process.in.text.trim()  // stndout
            return output.split("\\n")  // multiple rows may be returned
        } catch (Exception e) {
            throw new Exception("SQLite query failed:\nquery: $query\nException: $e\nCause: ${e.getCause()}")
        }
    }
}