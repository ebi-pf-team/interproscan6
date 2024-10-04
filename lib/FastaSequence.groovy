import java.security.MessageDigest

class FastaSequence implements Serializable {
    // DNA, RNA, gaps
    static final String NUCLEIC_ALPHABET = /(?i)[ACGTUN\-\.\s]+/
    // 20 standard AAs, Sec, Pyl, any/unknown, Asx, Glx, Xle, gaps
    static final String PROTEIN_ALPHABET = /(?i)[ACDEFGHIKLMNPQRSTVWYUOXBZJ\-\.\s]+/

    String id
    String description
    String sequence = ""
    String md5 = null
    FastaSequence translatedFrom = null

    FastaSequence(String id, String description) {
        this.id = id
        this.description = description
    }

    FastaSequence(String id, String description, String sequence) {
        this.id = id
        this.description = description
        this.sequence = sequence
        this.md5 = this.getMD5(sequence)
    }

    static String getMD5(String sequence) {
        MessageDigest md = MessageDigest.getInstance("MD5")
        byte[] digest = md.digest(sequence.toUpperCase().bytes)
        return digest.encodeHex().toString().toUpperCase()
    }

    void updateMD5() {
        this.md5 = this.getMD5(this.sequence)
    }

    static FastaSequence fromMap(Map data) {
        FastaSequence seq = new FastaSequence(data.id, data.description)
        seq.sequence = data.sequence
        seq.md5 = data.md5
        return seq
    }

    String validate(boolean isNucleic = false, String extraForbiddenChars = "") {
        def alphabet = isNucleic ? this.NUCLEIC_ALPHABET : this.PROTEIN_ALPHABET
        def pattern = "[^${alphabet}]"
        def invalidChars = this.sequence.findAll(/(?i)${pattern}/).unique().join("")
        if (invalidChars) {
            return invalidChars
        }

        if (!isNucleic && extraForbiddenChars) {
            pattern = /(?i)[\Q${extraForbiddenChars}\E]/
            invalidChars = this.sequence.findAll(pattern).unique().join("")
            if (invalidChars) {
                return invalidChars
            }
        }
        
        return null
    }
}
