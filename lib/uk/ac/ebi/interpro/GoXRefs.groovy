package uk.ac.ebi.interpro

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