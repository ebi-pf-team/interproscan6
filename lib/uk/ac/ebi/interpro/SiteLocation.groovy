package uk.ac.ebi.interpro

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