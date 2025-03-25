class SeqDatabase {
    // This class contains the methods for querying the internal IPS6 seq db, and handling any errors
    String dbPath

    SeqDatabase(String dbPath) {
        this.dbPath = dbPath
    }

    List<String> query(String query) {
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

    Map<String, List<String>> groupProteins(Map proteinMatches) {
        /* Gather nucleotide Seq IDs and child protein MD5s so we can gather all ORFs from the same
        parent NT seq together in the final output. */
        def nucleicRelationships = [:]  // [ntSeqId: [proteinMd5]]
        proteinMatches.each { proteinMd5, matchesNode ->
            // get all parent NT seq Ids
            def proteinQuery = """SELECT N.nt_md5
                FROM NUCLEOTIDE AS N
                LEFT JOIN PROTEIN_TO_NUCLEOTIDE AS N2P ON N.nt_md5 = N2P.nt_md5
                WHERE N2P.protein_md5 = '$proteinMd5'
                """
            def result = this.query(proteinQuery)
            result.each { ntSeqId ->
                nucleicRelationships.computeIfAbsent(ntSeqId, { [] as Set })
                nucleicRelationships[ntSeqId].add(proteinMd5)
            }
        }
        return nucleicRelationships
    }

    List<String> getSeqData(String md5, Boolean nucleic, String nucleicMd5 = null) {
        /* Retrieve all associated seq IDs and desc for the given seq (id'd by it's md5 hash)
        Return a list of rows as a seq may be associated with multiple ids */
        def seqQuery = ""
        if (nucleic && !nucleicMd5) { // retrieve data for a nucleic sequence
            seqQuery = """SELECT N.id, N.description, S.sequence
                   FROM NUCLEOTIDE AS N
                   LEFT JOIN NUCLEOTIDE_SEQUENCE AS S ON N.nt_md5 = S.nt_md5
                   WHERE N.nt_md5 = '$md5';"""
        } else if (nucleic && nucleicMd5) {  // retrieve the open reading frame for the current nucleic seq
            // a protein may be associated with multiple unique nucleotide sequences
            seqQuery = """SELECT P.id, P.description, S.sequence
                   FROM PROTEIN AS P
                   LEFT JOIN PROTEIN_SEQUENCE AS S ON P.protein_md5 = S.protein_md5
                   WHERE P.protein_md5 = '$nucleicMd5' AND P.description REGEXP 'source=${nucleicMd5}*'
                   """
        } else {  // retrieve data for a protein
            seqQuery = """SELECT P.id, P.description, S.sequence
                   FROM PROTEIN AS P
                   LEFT JOIN PROTEIN_SEQUENCE AS S ON P.protein_md5 = S.protein_md5
                   WHERE P.protein_md5 = '$md5';"""
        }
        return this.query(seqQuery)
    }
}