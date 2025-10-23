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