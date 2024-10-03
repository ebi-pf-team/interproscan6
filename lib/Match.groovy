class Match implements Serializable {
    String modelAccession
    Integer sequenceLength
    Double evalue
    Double score
    Double bias
    Signature signature = null
    List<Location> locations = []

    // PANTHER
    // String subfamilyAccession
    // String subfamilyName
    // String graftPoint
    // String proteinClass
    // String ancestralNodeID

    // SignalP
    // String orgType

    // PRINTS
    // String graphscan

    Match(String modelAccession) {
        this.modelAccession = modelAccession
    }

    Match(String modelAccession, Double evalue, Double score, Double bias) {
        this.modelAccession = modelAccession
        this.evalue = evalue
        this.score = score
        this.bias = bias
    }

    static Match fromMap(Map data) {
        Match match = new Match(data.modelAccession, data.evalue, data.score, data.bias)
        match.sequenceLength = data.sequenceLength
        match.signature = Signature.fromMap(data.signature)
        match.locations = data.locations.collect { Location.fromMap(it) }
        return match
    }

    void addLocation(Location location) {
        this.locations.add(location)
    }

    void setSequences(int locationIndex, String querySequence, String targetSequence) {
        Location location = this.locations[locationIndex]
        location.querySequence = querySequence
        location.targetSequence = targetSequence
    }
}

class Signature implements Serializable {
    String accession
    String name
    String description
    SignatureLibraryRelease signatureLibraryRelease
    Entry entry

    Signature(String accession, 
              String name, 
              String description, 
              SignatureLibraryRelease library,
              Entry entry) {
        this.accession = accession
        this.name = name
        this.description = description
        this.library = library
        this.entry = entry
    }

    static Signature fromMap(Map data) {
        return new Signature(
            data.accession,
            data.name,
            data.description,
            SignatureLibraryRelease.fromMap(data.signatureLibraryRelease),
            Entry.fromMap(data.entry)
        )
    }
}

class SignatureLibraryRelease implements Serializable {
    String library
    String version

    SignatureLibraryRelease(String library, String version) {
        this.library = library
        this.version = version
    }

    static SignatureLibraryRelease fromMap(Map data) {
        return new SignatureLibraryRelease(data.library, data.version)
    }
}

class Entry implements Serializable {
    String accession
    String name
    String description
    String type

    Entry(String accession, String name, String description, String type) {
        this.accession = accession
        this.name = name
        this.description = description
        this.type = type
    }

    static Entry fromMap(Map data) {
        return new Entry(data.accession, data.name, data.description, data.type)
    }
}

class Location implements Serializable {
    int start
    int end
    Integer hmmStart
    Integer hmmEnd
    Integer hmmLength
    String hmmBounds
    Integer envelopeStart
    Integer envelopeEnd
    Double evalue
    Double score
    Double bias
    String querySequence
    String targetSequence
    String sequenceFeature
    List<LocationFragment> fragments = []
    boolean representative = false

    // pvalue
    // level
    // cigarAlignment
    // sites
    // motifNumber

    Location(int start, 
             int end, 
             Integer hmmStart = null,
             Integer hmmEnd = null, 
             Integer hmmLength = null,
             String hmmBounds = null,
             Integer envelopeStart = null,
             Integer envelopeEnd = null,
             Double evalue = null,
             Double score = null,
             Double bias = null) {
        this.start = start
        this.end = end
        this.hmmStart = hmmStart
        this.hmmEnd = hmmEnd
        this.hmmLength = hmmLength
        this.hmmBounds = hmmBounds
        this.envelopeStart = envelopeStart
        this.envelopeEnd = envelopeEnd
        this.evalue = evalue
        this.score = score
        this.bias = bias
        LocationFragment fragment = new LocationFragment(start, end, "CONTINUOUS")
        this.fragments = [fragment]
    }

    Location(int start, int end, String sequenceFeature = null) {
        this.start = start
        this.end = end
        this.sequenceFeature = sequenceFeature
        LocationFragment fragment = new LocationFragment(start, end, "CONTINUOUS")
        this.fragments = [fragment]
    }

    static Location fromMap(data) {
        Location loc = new Location(
            data.start,
            data.end,
            data.hmmStart,
            data.hmmEnd,
            data.hmmLength,
            data.hmmBounds,
            data.envelopeStart,
            data.envelopeEnd,
            data.evalue,
            data.score,
            data.bias
        )
        loc.querySequence = data.querySequence
        loc.targetSequence = data.targetSequence
        loc.fragments = data.fragments.collect { LocationFragment.fromMap(it) }
        loc.representative = data.representative
        return loc
    }
}

class LocationFragment implements Serializable {
    int start
    int end
    String dcStatus

    LocationFragment(int start, int end, String dcStatus) {
        this.start = start
        this.end = end
        this.dcStatus = dcStatus
    }

    static LocationFragment fromMap(Map data) {
        return new LocationFragment(data.start, data.end, data.dcStatus)
    }
}
