import groovy.json.JsonOutput

import Match

process WRITE_FASTA {
    label    'tiny'
    executor 'local'

    input:
    tuple val(meta), val(fasta)

    output:
    tuple val(meta), path("sequences.fasta")

    exec:
    def outputFilePath = task.workDir.resolve("sequences.fasta")
    FastaFile.write(
        fasta.toString(),
        outputFilePath.toString(),
        "BJOZ",
        [:]
    )
}

process SEARCH_PHOBIUS {
    label       'small', 'ips6_container'
    stageInMode 'copy'

    input:
    tuple val(meta), path(fasta)
    path phobius_dir

    output:
    tuple val(meta), path("phobius.out")

    script:
    """
    sed -i 's|/decodeanhmm";|/decodeanhmm.64bit";|' ${phobius_dir}/phobius.pl
    ${phobius_dir}/phobius.pl ${fasta} > phobius.out
    rm -rf ${phobius_dir}
    """
}

process PARSE_PHOBIUS {
    label    'tiny'
    executor 'local'

    input:
    tuple val(meta), val(phobius_out)

    output:
    tuple val(meta), path("phobius.json")

    exec:
    def outputFilePath = task.workDir.resolve("phobius.json")
    def matches = [:]
    def tmpMatches = [:]
    def sequenceId = null
    boolean isSignalPeptide = false
    boolean isTransmembrane = false

    SignatureLibraryRelease library = new SignatureLibraryRelease("Phobius", "1.01")
    def signatures = [
        "CYTOPLASMIC_DOMAIN"     : new Signature("CYTOPLASMIC_DOMAIN", "Cytoplasmic domain", 
                                                 "Region of a membrane-bound protein predicted to be outside the membrane, in the cytoplasm", library, null),
        "NON_CYTOPLASMIC_DOMAIN" : new Signature("NON_CYTOPLASMIC_DOMAIN", "Non cytoplasmic domain", 
                                                 "Region of a membrane-bound protein predicted to be outside the membrane, in the extracellular region", library, null),
        "SIGNAL_PEPTIDE"         : new Signature("SIGNAL_PEPTIDE", "Signal Peptide", "Signal Peptide region", library, null),
        "SIGNAL_PEPTIDE_C_REGION": new Signature("SIGNAL_PEPTIDE_C_REGION", "Signal peptide C-region", 
                                                 "C-terminal region of a signal peptide", library, null),
        "SIGNAL_PEPTIDE_H_REGION": new Signature("SIGNAL_PEPTIDE_H_REGION", "Signal peptide H-region", 
                                                 "Hydrophobic region of a signal peptide", library, null),
        "SIGNAL_PEPTIDE_N_REGION": new Signature("SIGNAL_PEPTIDE_N_REGION", "Signal peptide N-region", 
                                                 "N-terminal region of a signal peptide", library, null),
        "TRANSMEMBRANE"          : new Signature("TRANSMEMBRANE", "Transmembrane region", 
                                                 "Region of a membrane-bound protein predicted to be embedded in the membrane", library, null),
    ]

    file(phobius_out.toString()).eachLine { line ->
        line = line.trim()
        if (line == "//") {
            assert sequenceId != null

            if (isSignalPeptide || isTransmembrane) {
                // Only consider sequences with predicted signal peptides or transmembrane helices
                matches[sequenceId] = tmpMatches.clone()  
            }

            sequenceId = null
            isSignalPeptide = false
            isTransmembrane = false
            tmpMatches.clear()
            return
        }
        def fields = line.trim().split(/\s+/, 5)
        String field = fields[0];
        if (field == "ID") {
            assert sequenceId == null
            assert fields.size() == 2
            sequenceId = fields[1]
        } else {
            assert field == "FT"
            assert fields.size() == 4 || fields.size() == 5
            String type = fields[1]
            int start = fields[2].toInteger()
            int end = fields[3].toInteger()
            String qualifier = null
            if (fields.size() == 5 && !fields[4].trim().isEmpty()) {
                qualifier = fields[4].trim()
            }

            String modelAccession = null
            if (type == "SIGNAL") {
                modelAccession = "SIGNAL_PEPTIDE"
                isSignalPeptide = true
            } else if (type == "DOMAIN") {
                switch (qualifier) {
                    case "CYTOPLASMIC.":
                        modelAccession = "CYTOPLASMIC_DOMAIN"
                        break
                    case "NON CYTOPLASMIC.":
                        modelAccession = "NON_CYTOPLASMIC_DOMAIN"
                        break
                    case "N-REGION.":
                        modelAccession = "SIGNAL_PEPTIDE_N_REGION"
                        isSignalPeptide = true
                        break
                    case "H-REGION.":
                        modelAccession = "SIGNAL_PEPTIDE_H_REGION"
                        isSignalPeptide = true
                        break
                    case "C-REGION.":
                        modelAccession = "SIGNAL_PEPTIDE_C_REGION"
                        isSignalPeptide = true
                        break
                }
            } else if (type == "TRANSMEM") {
                modelAccession = "TRANSMEMBRANE"
                isTransmembrane = true
            }
            
            if (modelAccession == null) {
                // Some features (e.g. REGION, TOPO_DOM) can be ignored
                return
            }

            Match match
            if (tmpMatches.containsKey(modelAccession)) {
                match = tmpMatches[modelAccession]
            } else {
                match = new Match(modelAccession)
                match.signature = signatures[modelAccession]
                tmpMatches[modelAccession] = match
            }

            match.addLocation(new Location(start, end))            
        }
    }

    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}