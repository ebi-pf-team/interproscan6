package uk.ac.ebi.interpro


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
            String residue = residueAnnotation.substring(0, 1)
            int position = residueAnnotation.substring(1).toInteger()
            siteLocations.add(new SiteLocation(residue, position, position))
        }
        return siteLocations
    }

    static Site fromMap(Map data) {
        Site site = new Site(
                data.description,
                data.siteLocations.collect { SiteLocation.fromMap(it) }
        )
    }

    boolean isInRange(int start, int end) {
        return start <= this.start && this.end <= end
    }

    @Override
    public int hashCode() {
        return Objects.hash(description, numLocations, siteLocations)
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true
        if (obj == null || getClass() != obj.getClass()) return false
        return (
                description == obj.description &&
                        numLocations == obj.numLocations &&
                        Objects.equals(siteLocations, obj.siteLocations)
        )
    }

    public Object clone() {
        return new Site(description, siteLocations.collect{ it.clone() })
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

        @Override
        public int hashCode() {
            return Objects.hash(residue, start, end)
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj) return true
            if (obj == null || getClass() != obj.getClass()) return false
            return (
                residue == obj.residue &&
                start == obj.start &&
                end == obj.end
            )
        }

        public Object clone() {
            return new SiteLocation(residue, start, end)
        }
    }
}