import groovy.json.JsonOutput

process SEARCH_PHOBIUS {
    label       'phobius_runner'
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
    label 'analysis_parser'

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

            if (!tmpMatches.containsKey(modelAccession)) {
                tmpMatches[modelAccession] = new Match(modelAccession)
            }

            tmpMatches[modelAccession].addLocation(new Location(start, end))            
        }
    }

    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}