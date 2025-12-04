package uk.ac.ebi.interpro

import uk.ac.ebi.interpro.GoXRefs
import uk.ac.ebi.interpro.PathwayXRefs

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
