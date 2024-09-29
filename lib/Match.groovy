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

    static class Signature implements Serializable {
        String accession
        String name
        String description
        SignatureLibraryRelease signatureLibraryRelease
        Entry entry
    }

    static class SignatureLibraryRelease implements Serializable {
        String library
        String version
    }

    static class Entry implements Serializable {
        String accession
        String name
        String description
        String type
    }

    static class Location implements Serializable {
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
        String alignment
        String sequenceFeature
        List<LocationFragment> fragments
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
                 Double bias = null,
                 String alignment = null,
                 String sequenceFeature = null,
                 List<LocationFragment> fragments = null) {
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
            this.alignment = alignment
            this.sequenceFeature = sequenceFeature
            if (fragments) {
                this.fragments = fragments
            } else {
                LocationFragment fragment = new LocationFragment(start, end, "CONTINUOUS")
                this.fragments = [fragment]
            }
        }

        Location(int start, int end, String sequenceFeature) {
            this.start = start
            this.end = end
            this.sequenceFeature = sequenceFeature
            LocationFragment fragment = new LocationFragment(start, end, "CONTINUOUS")
            this.fragments = [fragment]
        }
    }

    static class LocationFragment implements Serializable {
        int start
        int end
        String dcStatus

        LocationFragment(int start, int end, String dcStatus) {
            this.start = start
            this.end = end
            this.dcStatus = dcStatus
        }
    }

    Match(String modelAccession) {
        this.modelAccession = modelAccession
    }

    Match(String modelAccession, Double evalue, Double score, Double bias) {
        this.modelAccession = modelAccession
        this.evalue = evalue
        this.score = score
        this.bias = bias
    }

    void addMobiDBLiteLocation(int start, int end, String sequenceFeature) {
        Location location = new Location(start, end, sequenceFeature)
        this.locations.add(location)
    }
}