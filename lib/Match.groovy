class Match implements Serializable {
    String modelAccession
    Integer sequenceLength
    Double evalue
    Double score
    Double bias
    Signature signature = null
    List<Location> locations = []
    boolean included = true  // for HMMER3 matches (inclusion threshold)
    RepresentativeInfo representativeInfo = null

    // PANTHER
    TreeGrafter treegrafter = null

    // SignalP
    SignalP signalp = null

    // PRINTS
    // String graphscan

    Match(String modelAccession) {
        this.modelAccession = modelAccession
    }

    Match(String modelAccession, Double evalue, Double score) {
        this.modelAccession = modelAccession
        this.evalue = evalue
        this.score = score
    }

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
        match.included = data.included
        match.locations = data.locations.collect { Location.fromMap(it) }
        match.treegrafter = TreeGrafter.fromMap(data.treegrafter)
        return match
    }

    void addLocation(Location location) {
        this.locations.add(location)
    }

    void addSignalPeptide(String orgType, int cleavageSiteStart, int cleavageSiteEnd) {
        this.signalp = new SignalP(orgType, cleavageSiteStart, cleavageSiteEnd)
    }

    void setAlignments(int locationIndex, String queryAlignment, String targetAlignment) {
        Location location = this.locations[locationIndex]
        location.queryAlignment = queryAlignment
        location.targetAlignment = targetAlignment
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
    List<GoXRefs> goXRefs = []
    List<PathwayXRefs> pathwayXRefs = []

    Entry(String accession,
          String name,
          String description,
          String type) {
        this.accession = accession
        this.name = name
        this.description = description
        this.type = type
    }

    Entry(String accession,
          String name,
          String description,
          String type,
          List<GoXRefs> goXRefs,
          List<PathwayXRefs> pathwayXRefs) {
        this.accession = accession
        this.name = name
        this.description = description
        this.type = type
        this.goXRefs = goXRefs
        this.pathwayXRefs = pathwayXRefs
    }

    static Entry fromMap(Map data) {
        if (data == null) {
            return null
        }
        return new Entry(
                data.accession,
                data.name,
                data.description,
                data.type,
                data.goXRefs.collect { GoXRefs.fromMap(it) },
                data.pathwayXRefs.collect { PathwayXRefs.fromMap(it) }
        )
    }

    void addGoXRefs(GoXRefs go) {
        this.goXRefs.add(go)
    }

    void addPathwayXRefs(PathwayXRefs pa) {
        this.pathwayXRefs.add(pa)
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
    String queryAlignment
    String targetAlignment
    String sequenceFeature
    String level
    String cigarAlignment
    List<LocationFragment> fragments = []
    List<Site> sites = []
    boolean representative = false
    boolean included = true  // for HMMER3 matches (inclusion threshold)
    Float pvalue = null // SignalP

    // pvalue
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

    Location(int start, int end, Double score, String targetAlignment) {
        this.start = start
        this.end = end
        this.score = score
        LocationFragment fragment = new LocationFragment(start, end, "CONTINUOUS")
        this.fragments = [fragment]
        this.targetAlignment = targetAlignment
    }

    Location(int start, int end, Double evalue, List<LocationFragment> fragments) {
        this.start = start
        this.end = end
        this.evalue = evalue
        this.fragments = fragments
    }

    Location(int start, int end, String level, String alignment, String cigarAlignment) {
        this.start = start
        this.end = end
        this.level = level
        LocationFragment fragment = new LocationFragment(start, end, "CONTINUOUS")
        this.fragments = [fragment]
        this.targetAlignment = alignment
        this.cigarAlignment = cigarAlignment
    }

    Location(int start, int end, float pvalue) { // Used for SignalP
        this.start = start
        this.end = end
        this.pvalue = pvalue
        LocationFragment fragment = new LocationFragment(start, end, "CONTINUOUS")
        this.fragments = [fragment]
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
        loc.queryAlignment = data.queryAlignment
        loc.targetAlignment = data.targetAlignment
        loc.fragments = data.fragments.collect { LocationFragment.fromMap(it) }
        loc.representative = data.representative
        loc.included = data.included
        loc.sites = data.sites
        loc.sequenceFeature = data.sequenceFeature
        loc.level = data.level
        loc.cigarAlignment = data.cigarAlignment
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
    String label = null
    String group = null
    int hmmStart
    int hmmEnd
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

    // PIRSR case
    Site(String description,
        int group,
        int hmmEnd,
        int hmmStart,
        String label,
        List<SiteLocation> siteLocations) {
        this.description = description
        this.group = group
        this.hmmEnd = hmmEnd
        this.hmmStart = hmmStart
        this.label = label
        this.numLocations = siteLocations.size()
        this.siteLocations = siteLocations
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

class TreeGrafter implements Serializable {
    String ancestralNodeID
    String graftPoint
    String subfamilyAccession
    String subfamilyName
    String subfamilyDescription
    String proteinClass

    TreeGrafter(String ancestralNodeID) {
        this.ancestralNodeID = ancestralNodeID
    }

    static TreeGrafter fromMap(Map data) {
        if (data == null) {
            return null
        }
        TreeGrafter tg = new TreeGrafter(data.ancestralNodeID)
        tg.graftPoint = data.graftPoint
        tg.subfamilyAccession = data.subfamilyAccession
        tg.subfamilyName = data.subfamilyName
        tg.subfamilyDescription = data.subfamilyDescription
        tg.proteinClass = data.proteinClass
        return tg
    }
}

class SignalP implements Serializable {
    String orgType
    int cleavageSiteStart
    int cleavageSiteEnd

    SignalP(String orgType, int cleavageSiteStart, int cleavageSiteEnd) {
        this.orgType = orgType
        this.cleavageSiteStart = cleavageSiteStart
        this.cleavageSiteEnd = cleavageSiteEnd
    }

    static SignalP fromMap(Map data) {
        if (data == null) {
            return null
        }
        SignalP sp = new SignalP(
                data.orgType,
                data.cleavageSiteStart,
                data.cleavageSiteEnd
        )
        return sp
    }
}

class RepresentativeInfo implements Serializable {
    String type
    int rank

    RepresentativeInfo(String type, int rank) {
        this.type = type
        this.rank = rank
    }

    static RepresentativeInfo fromMap(Map data) {
        if (data == null) {
            return null
        }
        return new RepresentativeInfo(data.type, data.rank)
    }
}

class GoXRefs implements Serializable {
    String name
    String databaseName
    String category
    String id

    GoXRefs(String name, String databaseName, String category, String id) {
        this.name = name
        this.databaseName = databaseName
        this.category = category
        this.id = id
    }

    static GoXRefs fromMap(Map data) {
        if (data == null) {
            return null
        }
        return new GoXRefs(data.name, data.databaseName, data.category, data.id)
    }
}

class PathwayXRefs implements Serializable {
    String name
    String databaseName
    String id

    PathwayXRefs(String name, String databaseName, String id) {
        this.name = name
        this.databaseName = databaseName
        this.id = id
    }

    static PathwayXRefs fromMap(Map data) {
        if (data == null) {
            return null
        }
        return new PathwayXRefs(data.name, data.databaseName, data.id)
    }
}
