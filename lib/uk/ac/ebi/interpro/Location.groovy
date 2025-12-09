package uk.ac.ebi.interpro
import uk.ac.ebi.interpro.LocationFragment

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
    Double pvalue
    Integer motifNumber
    String level
    String cigarAlignment
    List<LocationFragment> fragments = []
    List<Site> sites = []
    boolean representative = false
    boolean included = true  // for HMMER3 matches (inclusion threshold)

    Location() {}

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

    Location(int start,
             int end,
             Integer hmmStart,
             Integer hmmEnd,
             Integer envelopeStart,
             Integer envelopeEnd,
             Double evalue,
             Double score,
             Double bias,
             List<LocationFragment> fragments) { // Used for SFLD
        this.start = start
        this.end = end
        this.hmmStart = hmmStart
        this.hmmEnd = hmmEnd
        this.envelopeStart = envelopeStart
        this.envelopeEnd = envelopeEnd
        this.evalue = evalue
        this.score = score
        this.bias = bias
        this.fragments = fragments
    }

    Location(int start, int end, String sequenceFeature = null) { // Used for Coils, MobiDB, Phobius
        this.start = start
        this.end = end
        this.sequenceFeature = sequenceFeature
        LocationFragment fragment = new LocationFragment(start, end, "CONTINUOUS")
        this.fragments = [fragment]
    }

    Location(int start, int end, Double evalue, Double score) { // Used for CDD
        this.start = start
        this.end = end
        this.evalue = evalue
        this.score = score
        LocationFragment fragment = new LocationFragment(start, end, "CONTINUOUS")
        this.fragments = [fragment]
    }

    Location(int start, int end, Double score, String targetAlignment, cigarAlignment) { // Used for Hamap, PrositeProfiles
        this.start = start
        this.end = end
        this.score = score
        LocationFragment fragment = new LocationFragment(start, end, "CONTINUOUS")
        this.fragments = [fragment]
        this.targetAlignment = targetAlignment
        this.cigarAlignment = cigarAlignment
    }

    Location(int start, int end, Double pvalue, Double score, Integer motifNumber) { // Used for PRINTS
        this.start = start
        this.end = end
        this.pvalue = pvalue
        this.score = score
        this.motifNumber = motifNumber
        LocationFragment fragment = new LocationFragment(start, end, "CONTINUOUS")
        this.fragments = [fragment]
    }

    Location(int start, int end, Integer hmmLength, Double evalue, List<LocationFragment> fragments) { // Used for Superfamily
        this.start = start
        this.end = end
        this.hmmLength = hmmLength
        this.evalue = evalue
        this.fragments = fragments
    }

    Location(int start, int end, Double score, List<LocationFragment> fragments) { // Used for InterPro-N
        this.start = start
        this.end = end
        this.score = score
        this.fragments = fragments
    }

    Location(int start, int end, String level, String targetAlignment, String cigarAlignment) { // Used for PrositePatterns
        this.start = start
        this.end = end
        this.level = level
        LocationFragment fragment = new LocationFragment(start, end, "CONTINUOUS")
        this.fragments = [fragment]
        this.targetAlignment = targetAlignment
        this.cigarAlignment = cigarAlignment
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
        loc.pvalue = data.pvalue
        loc.motifNumber = data.motifNumber
        return loc
    }

    static String getHmmBounds(String hmmBounds) {
        return [
                "[]"  : "COMPLETE",
                "[."  : "N_TERMINAL_COMPLETE",
                ".]"  : "C_TERMINAL_COMPLETE",
                ".."  : "INCOMPLETE"
        ][hmmBounds]
    }

    static String getReverseHmmBounds(String hmmBounds) {
        return [
                "COMPLETE"            : "[]",
                "N_TERMINAL_COMPLETE" : "[.",
                "C_TERMINAL_COMPLETE" : ".]",
                "INCOMPLETE"          : ".."
        ][hmmBounds]
    }

    @Override
    public int hashCode() {
        return Objects.hash(start, end, hmmStart, hmmEnd, hmmLength, hmmBounds,
            envelopeStart, envelopeEnd, evalue, score, bias,
            queryAlignment, targetAlignment, fragments, sites)
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true
        if (obj == null || getClass() != obj.getClass()) return false
        return (
            start == obj.start &&
            end == obj.end &&
            hmmStart == obj.hmmStart &&
            hmmEnd == obj.hmmEnd &&
            hmmLength == obj.hmmLength &&
            hmmBounds == obj.hmmBounds &&
            envelopeStart == obj.envelopeStart &&
            envelopeEnd == obj.envelopeEnd &&
            evalue == obj.evalue &&
            score == obj.score &&
            bias == obj.bias &&
            queryAlignment == obj.queryAlignment &&
            targetAlignment == obj.targetAlignment &&
            Objects.equals(fragments, obj.fragments) &&
            Objects.equals(sites, obj.sites)
        )
    }

    public Object clone() {
        Location loc = new Location(
            start,
            end,
            hmmStart,
            hmmEnd,
            hmmLength,
            hmmBounds,
            envelopeStart,
            envelopeEnd,
            evalue,
            score,
            bias
        )
        loc.queryAlignment = queryAlignment
        loc.targetAlignment = targetAlignment
        loc.fragments = fragments.collect { it.clone() }
        loc.representative = representative
        loc.included = included
        loc.sites = sites.collect { it.clone() }
        return loc
    }
}