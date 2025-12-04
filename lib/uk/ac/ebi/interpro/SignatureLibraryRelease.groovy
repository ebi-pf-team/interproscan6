package uk.ac.ebi.interpro

class SignatureLibraryRelease implements Serializable {
    String library
    String version

    SignatureLibraryRelease(String library, String version) {
        this.library = library
        this.version = version
    }

    static SignatureLibraryRelease fromMap(Map data) {
        if (data == null) {
            return null
        }
        return new SignatureLibraryRelease(data.library, data.version)
    }
}