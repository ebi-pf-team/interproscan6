package uk.ac.ebi.interpro

import uk.ac.ebi.interpro.Location
import uk.ac.ebi.interpro.RepresentativeInfo
import uk.ac.ebi.interpro.Signature
import uk.ac.ebi.interpro.TreeGrafter

class Match implements Serializable {
    String modelAccession
    Double evalue
    Double score
    Double bias
    Signature signature = null
    List<Location> locations = []
    boolean included = true  // for HMMER3 matches (inclusion threshold)
    RepresentativeInfo representativeInfo = null
    String source = null

    // PANTHER
    TreeGrafter treegrafter = null

    // PRINTS
    String graphscan = null

    Match(String modelAccession) {
        this.modelAccession = modelAccession
    }

    Match(String modelAccession, Signature signature) {
        this.modelAccession = modelAccession
        this.signature = signature
        this.source = signature.signatureLibraryRelease.library
    }

    Match(String modelAccession, Double evalue, String graphscan, Signature signature) {
        this.modelAccession = modelAccession
        this.evalue = evalue
        this.graphscan = graphscan
        this.signature = signature
        this.source = signature.signatureLibraryRelease.library
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

    Match(String modelAccession, Double evalue, Double score, Signature signature) {
        this.modelAccession = modelAccession
        this.evalue = evalue
        this.score = score
        this.signature = signature
        this.source = signature.signatureLibraryRelease.library
    }

    Match(String modelAccession, Double evalue, Double score, Double bias, Signature signature) {
        this.modelAccession = modelAccession
        this.evalue = evalue
        this.score = score
        this.bias = bias
        this.signature = signature
        this.source = signature.signatureLibraryRelease.library
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
        Match match = new Match(data.modelAccession)
        match.evalue = data.evalue
        match.score = data.score
        match.bias = data.bias
        match.signature = Signature.fromMap(data.signature)
        match.included = data.included
        match.locations = data.locations.collect { Location.fromMap(it) }
        match.treegrafter = TreeGrafter.fromMap(data.treegrafter)
        match.graphscan = data.graphscan
        match.representativeInfo = RepresentativeInfo.fromMap(data.representativeInfo)
        match.source = data.source
        return match
    }

    void addLocation(Location location) {
        this.locations.add(location)
    }

    void setAlignments(int locationIndex, String queryAlignment, String targetAlignment) {
        Location location = this.locations[locationIndex]
        location.queryAlignment = queryAlignment
        location.targetAlignment = targetAlignment
    }

    static String encodeCigarAlignment(String alignment) {
        if (!alignment) return ""
        StringBuilder cigar = new StringBuilder()
        char prevOp = '\0'
        int count = 0

        for (int i = 0; i < alignment.length(); i++) {
            char c = alignment.charAt(i)
            char op

            if (c == '-') {
                op = 'D'
            } else if (Character.isLowerCase(c)) {
                op = 'I'
            } else {
                op = 'M'
            }

            if (op == prevOp) {
                count++
            } else {
                if (count > 0) {
                    cigar.append(count).append(prevOp)
                }
                prevOp = op
                count = 1
            }
        }

        // Append the final operation
        if (count > 0) {
            cigar.append(count).append(prevOp)
        }

        return cigar.toString()
    }

    @Override
    public int hashCode() {
        int x = Objects.hash(modelAccession, evalue, score, bias, signature, locations)
        return x
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true
        if (obj == null || getClass() != obj.getClass()) return false
        return (
            modelAccession == obj.modelAccession &&
            evalue == obj.evalue &&
            score == obj.score &&
            bias == obj.bias &&
            signature == obj.signature &&
            locations == obj.locations
        )
    }
}
