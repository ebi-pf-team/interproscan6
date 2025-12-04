package uk.ac.ebi.interpro

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