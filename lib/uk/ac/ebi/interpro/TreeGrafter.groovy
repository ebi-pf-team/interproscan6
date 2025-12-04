package uk.ac.ebi.interpro

import uk.ac.ebi.interpro.GoXRefs

class TreeGrafter implements Serializable {
    String ancestralNodeID
    String graftPoint
    String subfamilyAccession
    String subfamilyName
    String subfamilyDescription
    String proteinClass
    List<GoXRefs> goXRefs = []

    TreeGrafter(String ancestralNodeID) {
        this.ancestralNodeID = ancestralNodeID
    }

    static TreeGrafter fromMap(Map data) {
        if (data == null) {
            return null
        }
        TreeGrafter tg = new TreeGrafter(data.ancestralNodeID)
        tg.graftPoint = data.graftPoint
        tg.subfamilyAccession = data.subfamilyAccession
        tg.subfamilyName = data.subfamilyName
        tg.subfamilyDescription = data.subfamilyDescription
        tg.proteinClass = data.proteinClass
        tg.goXRefs = data.goXRefs.collect { GoXRefs.fromMap(it) }
        return tg
    }

    void addGoXRefs(GoXRefs go) {
        this.goXRefs.add(go)
    }
}