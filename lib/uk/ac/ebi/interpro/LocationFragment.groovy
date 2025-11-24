package uk.ac.ebi.interpro

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

    @Override
    public int hashCode() {
        return Objects.hash(start, end, dcStatus)
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true
        if (obj == null || getClass() != obj.getClass()) return false
        return (
            start == obj.start &&
            end == obj.end &&
            dcStatus == obj.dcStatus
        )
    }

    public Object clone() {
        return new LocationFragment(start, end, dcStatus)
    }
}