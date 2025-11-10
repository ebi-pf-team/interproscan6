package uk.ac.ebi.interpro

class Signature implements Serializable {
    // The order of the fields here determines their order in the final output files
    String accession
    String name
    String description
    String type
    SignatureLibraryRelease signatureLibraryRelease = new SignatureLibraryRelease(null, null)
    Entry entry = null

    Signature(String accession) {
        this.accession = accession
    }

    Signature(String accession, SignatureLibraryRelease library) {
        this.accession = accession
        this.signatureLibraryRelease = library
    }

    Signature(String accession,
              String name,
              String description,
              SignatureLibraryRelease library,
              Entry entry) {
        this.accession = accession
        this.name = name
        this.description = description
        this.signatureLibraryRelease = library
        this.entry = entry
    }

    Signature(String accession,
              String name,
              String description,
              String type,
              SignatureLibraryRelease library,
              Entry entry) {
        this.accession = accession
        this.name = name
        this.description = description
        this.type = type
        this.signatureLibraryRelease = library
        this.entry = entry
    }

    void setType(String type) {
        this.type = type
    }

    static Signature fromMap(Map data) {
        if (data == null) {
            return null
        }
        return new Signature(
                data.accession,
                data.name,
                data.description,
                data.containsKey("type") ? data.type : null,  // Provide a default value (null) if 'type' is missing
                SignatureLibraryRelease.fromMap(data.signatureLibraryRelease),
                Entry.fromMap(data.entry)
        )
    }
}