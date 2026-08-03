process WRITE_FASTA {
    label    'mem_min'
    label    'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(fasta)

    output:
    tuple val(meta), path("sequences.fasta")

    exec:
    uk.ac.ebi.interpro.FastaFile.write(
        fasta,
        task.workDir.resolve("sequences.fasta"),
        "BJOZ",
        [:]
    )
}

process SEARCH_PHOBIUS {
    label       'mem_min'
    label       'time_short'
    stageInMode 'copy'
    container   'interpro/phobius:1.01'

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
    label    'mem_low'
    label    'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(phobius_out)

    output:
    tuple val(meta), path("phobius.json")

    exec:
    def matches = [:]
    def tmpMatches = [:]
    def sequenceId = null
    def isSignalPeptide = false
    def isTransmembrane = false
    def library = new uk.ac.ebi.interpro.SignatureLibraryRelease("Phobius", "1.01")
    def signatures = [
        "CYTOPLASMIC_DOMAIN"     : new uk.ac.ebi.interpro.Signature("CYTOPLASMIC_DOMAIN", "Cytoplasmic domain", 
                                                 "Region of a membrane-bound protein predicted to be outside the membrane, in the cytoplasm", "Region", library, null),
        "NON_CYTOPLASMIC_DOMAIN" : new uk.ac.ebi.interpro.Signature("NON_CYTOPLASMIC_DOMAIN", "Non cytoplasmic domain", 
                                                 "Region of a membrane-bound protein predicted to be outside the membrane, in the extracellular region", "Region", library, null),
        "SIGNAL_PEPTIDE"         : new uk.ac.ebi.interpro.Signature("SIGNAL_PEPTIDE", "Signal Peptide", "Signal Peptide region", "Region", library, null),
        "SIGNAL_PEPTIDE_C_REGION": new uk.ac.ebi.interpro.Signature("SIGNAL_PEPTIDE_C_REGION", "Signal peptide C-region", 
                                                 "C-terminal region of a signal peptide", "Region", library, null),
        "SIGNAL_PEPTIDE_H_REGION": new uk.ac.ebi.interpro.Signature("SIGNAL_PEPTIDE_H_REGION", "Signal peptide H-region", 
                                                 "Hydrophobic region of a signal peptide", "Region", library, null),
        "SIGNAL_PEPTIDE_N_REGION": new uk.ac.ebi.interpro.Signature("SIGNAL_PEPTIDE_N_REGION", "Signal peptide N-region", 
                                                 "N-terminal region of a signal peptide", "Region", library, null),
        "TRANSMEMBRANE"          : new uk.ac.ebi.interpro.Signature("TRANSMEMBRANE", "Transmembrane region", 
                                                 "Region of a membrane-bound protein predicted to be embedded in the membrane", "Region", library, null),
    ]

    phobius_out.eachLine { line ->
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
        def field = fields[0];
        if (field == "ID") {
            assert sequenceId == null
            assert fields.size() == 2
            sequenceId = fields[1]
        } else {
            assert field == "FT"
            assert fields.size() == 4 || fields.size() == 5
            def type = fields[1]
            def start = fields[2].toInteger()
            def end = fields[3].toInteger()
            def qualifier = null
            if (fields.size() == 5 && !fields[4].trim().isEmpty()) {
                qualifier = fields[4].trim()
            }

            def modelAccession = null
                        if (type == "SIGNAL") {
                modelAccession = "SIGNAL_PEPTIDE"
                isSignalPeptide = true
            } else if (type == "DOMAIN") {
                if (qualifier == "CYTOPLASMIC.") {
                    modelAccession = "CYTOPLASMIC_DOMAIN"
                } else if (qualifier == "NON CYTOPLASMIC.") {
                    modelAccession = "NON_CYTOPLASMIC_DOMAIN"
                } else if (qualifier == "N-REGION.") {
                    modelAccession = "SIGNAL_PEPTIDE_N_REGION"
                    isSignalPeptide = true
                } else if (qualifier == "H-REGION.") {
                    modelAccession = "SIGNAL_PEPTIDE_H_REGION"
                    isSignalPeptide = true
                } else if (qualifier == "C-REGION.") {
                    modelAccession = "SIGNAL_PEPTIDE_C_REGION"
                    isSignalPeptide = true
                }
            } else if (type == "TRANSMEM") {
                modelAccession = "TRANSMEMBRANE"
                isTransmembrane = true
            }
            
            if (modelAccession == null) {
                // Some features (e.g. REGION, TOPO_DOM) can be ignored
                return
            }

            def match
            if (tmpMatches.containsKey(modelAccession)) {
                match = tmpMatches[modelAccession]
            } else {
                def signature = signatures[modelAccession]
                match = new uk.ac.ebi.interpro.Match(modelAccession, signature)
                tmpMatches[modelAccession] = match
            }

            match.addLocation(new uk.ac.ebi.interpro.Location(start, end))            
        }
    }

    def filepath = task.workDir.resolve("phobius.json")
    filepath.text  = groovy.json.JsonOutput.toJson(matches)
}