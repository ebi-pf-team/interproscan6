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

    // Used for HAMAP, MobiDB-lite
    Match(String modelAccession) {
        this.modelAccession = modelAccession
    }

    // Used for CDD
    Match(String modelAccession, Double evalue, Double score) {
        this.modelAccession = modelAccession
        this.evalue = evalue
        this.score = score
    }

    // Used for Coils
    Match(String modelAccession, Double evalue, Double score, Double bias) {
        this.modelAccession = modelAccession
        this.evalue = evalue
        this.score = score
        this.bias = bias
    }

    void addSite(Site site) {
        for (Location location: this.locations) {
            if (site.isInRange(location.start, location.end)) {
                location.addSite(site)
                return
            }
        }
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

    Signature(String accession) {
        this.accession = accession
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

    static Signature fromMap(Map data) {
        if (data == null) {
            return null
        }
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
        if (data == null) {
            return null
        }
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
        if (data == null) {
            return null
        }
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
    List<Site> sites = []
    boolean representative = false

    // pvalue
    // level
    // cigarAlignment
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

    Location(int start, 
             int end, 
             Integer hmmStart,
             Integer hmmEnd, 
             Integer hmmLength,
             String hmmBounds,
             Integer envelopeStart,
             Integer envelopeEnd,
             Double evalue,
             Double score,
             Double bias,
             List<LocationFragment> fragments) {
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
        this.fragments = fragments
    }

    Location(int start, int end, String sequenceFeature = null) {
        this.start = start
        this.end = end
        this.sequenceFeature = sequenceFeature
        LocationFragment fragment = new LocationFragment(start, end, "CONTINUOUS")
        this.fragments = [fragment]
    }

    Location(int start, int end, Double score, String alignment) {
        this.start = start
        this.end = end
        this.score = score
        LocationFragment fragment = new LocationFragment(start, end, "CONTINUOUS")
        this.fragments = [fragment]
        this.targetSequence = alignment
    }

    void addSite(Site site) {
        this.sites.add(site)
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
        loc.sites = data.sites
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

class Site implements Serializable {
    String description
    int numLocations
    List<SiteLocation> siteLocations = []
    private int start = -1
    private int end = -1

    Site(String description, List<SiteLocation> siteLocations) {
        this.description = description
        this.siteLocations = siteLocations
        this.numLocations = siteLocations.size()

        for (SiteLocation loc: siteLocations) {
            if (this.start == -1 || loc.start < this.start) {
                this.start = loc.start
            }

            if (this.end == -1 || loc.end > this.end) {
                this.end = loc.end
            }
        }
    }

    Site(String description, String residues) {
        this(description, Site.getSiteLocationsFromString(residues))        
    }

    private static List<SiteLocation> getSiteLocationsFromString(String residues) {
        def residueAnnotations = residues.split(",")
        List<SiteLocation> siteLocations = []
        for (String residueAnnotation: residueAnnotations) {
            String residue = residueAnnotation.substring(0, 1);
            int position = residueAnnotation.substring(1).toInteger()
            siteLocations.add(new SiteLocation(residue, position, position))
        }
        return siteLocations
    }

    static Site fromMap(Map data) {
        return new Site(
            data.description,
            data.siteLocations.collect { SiteLocation.fromMap(it) })
    }

    boolean isInRange(int start, int end) {
        return start <= this.start && this.end <= end
    }
}

class SiteLocation implements Serializable {
    String residue
    int start
    int end

    SiteLocation(String residue, int start, int end) {
        this.start = start
        this.end = end
        this.residue = residue
    }

    static SiteLocation fromMap(Map data) {
        return new SiteLocation(data.start, data.end, data.residue)
    }
}
