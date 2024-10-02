import java.security.MessageDigest

class FastaSequence implements Serializable {
    // DNA, RNA, gaps
    static final String NUCLEIC_ALPHABET = ~/(?i)[ACGTUN\-\.\s]+/
    // 20 standard AAs, Sec, Pyl, any/unknown, Asx, Glx, Xle, gaps
    static final String PROTEIN_ALPHABET = /(?i)[ACDEFGHIKLMNPQRSTVWYUOXBZJ\-\.\s]+/

    String id
    String description
    String sequence
    String md5 = null

    void setMD5() {
        MessageDigest md = MessageDigest.getInstance("MD5")
        byte[] digest = md.digest(sequence.toUpperCase().bytes)
        this.md5 = digest.encodeHex().toString().toUpperCase()
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
