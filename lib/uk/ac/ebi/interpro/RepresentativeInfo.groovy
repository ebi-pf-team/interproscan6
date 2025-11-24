package uk.ac.ebi.interpro

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