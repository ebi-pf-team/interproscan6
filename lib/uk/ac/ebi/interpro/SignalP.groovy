package uk.ac.ebi.interpro

class SignalP {
    // Library reported for both organism variants
    static final String LIBRARY = "SignalP"

    // The only mode the Matches API stores pre-calculated results for
    static final String API_MODE = "fast"

    // SignalP --organism option -> appsConfig label
    static final Map<String, String> ORGANISMS = [
        eukarya: "signalp_euk",
        other  : "signalp_prok",
    ]

    static final List<String> APPLICATIONS = ORGANISMS.values() as List

    // Model accessions encode the mode and organism, e.g. "SignalP_fast_eukarya"
    static String modelAccession(String mode, String organism) {
        return "SignalP_${mode}_${organism}".toString()
    }

    static String mode(String modelAccession) {
        def parts = modelAccession.split("_")
        assert parts.size() == 3
        return parts[1]
    }

    static String organism(String modelAccession) {
        def parts = modelAccession.split("_")
        assert parts.size() == 3
        return parts[2]
    }

    // "SignalP_fast_eukarya" -> "signalp_euk"; null when unrecognised
    static String toLabel(String modelAccession) {
        return ORGANISMS[organism(modelAccession)]
    }
}
