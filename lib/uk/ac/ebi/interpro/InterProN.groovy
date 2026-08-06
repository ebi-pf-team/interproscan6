package uk.ac.ebi.interpro

class InterProN {
    // appsConfig label -> library name as reported by InterPro
    static final Map<String, String> SUPPORTED_DATABASES = [
        cathgene3d     : "CATH-Gene3D",
        cdd            : "CDD",
        hamap          : "HAMAP",
        ncbifam        : "NCBIFAM",
        panther        : "PANTHER",
        pfam           : "Pfam",
        pirsf          : "PIRSF",
        prints         : "PRINTS",
        prositeprofiles: "PROSITE profiles",
        prositepatterns: "PROSITE patterns",
        sfld           : "SFLD",
        smart          : "SMART",
        superfamily    : "SUPERFAMILY",
    ]

    // Library names emitted by the InterPro-N model that differ from the appsConfig label
    static final Map<String, String> LIBRARY_ALIASES = [ssf: "superfamily"]

    // Normalise a library or application name to an appsConfig label
    static String toLabel(String name) {
        def key = InterProScan.normalizeName(name)
        return LIBRARY_ALIASES.get(key, key)
    }

    // Supported databases to report, given the user's requested applications
    static Set<String> selectDatabases(List<String> applications) {
        def overlap = SUPPORTED_DATABASES.keySet().intersect(applications as Set)
        return (overlap ?: SUPPORTED_DATABASES.keySet()) as Set
    }
}
