process SEARCH_CDD {
    label     'mem_min', 'time_short'
    container 'interpro/cdd:1.0'

    input:
    tuple val(meta), path(fasta)
    tuple path(dir), val(rpsblast_db), val(rpsproc_db)

    output:
    tuple val(meta), path("rpsbproc.out")

    script:
    """
    # Find rpsblast location and set appropriate lib path
    RPSBLAST_PATH=\$(which rpsblast)
    BLAST_DIR=\$(dirname "\$RPSBLAST_PATH")
    if [ -d "\${BLAST_DIR}/lib" ]; then
        export LD_LIBRARY_PATH="\${BLAST_DIR}/lib:\$LD_LIBRARY_PATH"
    fi

    rpsblast \
        -query ${fasta} \
        -db "${dir}/${rpsblast_db}" \
        -out rpsblast.out \
        -evalue 0.01 -seg no -outfmt 11

    rpsbproc \
        --infile rpsblast.out \
        --outfile rpsbproc.out \
        --data-path ${dir}/${rpsproc_db} \
        -m std
    """
}

process PARSE_CDD {
    label    'mem_low', 'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(rpsbproc_out)

    output:
    tuple val(meta), path("cdd.json")

    exec:
    def library = new uk.ac.ebi.interpro.SignatureLibraryRelease("CDD", null)
    def sessionId = null
    def sequenceId = null
    def inDomains = false
    def inSites = false
    def hits = [:]
    def pssmHits = [:]
    def pssmSites = [:]
    rpsbproc_out.eachLine { line -> 
        if (line.startsWith("SESSION")) {
            // #SESSION        <session-ordinal>       <program>       <database>      <score-matrix>  <evalue-threshold>
            sessionId = line.split("\t")[1]
        } else if (line.startsWith("QUERY")) {
            // #QUERY  <query-id>      <seq-type>      <seq-length>    <definition-line>
            def fields = line.split("\t", 5)
            def definition = fields[4]
            // Only keep the sequence identifier
            sequenceId = definition.split()[0]
            hits[sequenceId] = [:]
        } else if (line.startsWith("DOMAINS")) {
            assert inDomains == false
            assert inSites == false
            inDomains = true
        } else if (line.startsWith("SITES")) {
            assert inDomains == false
            assert inSites == false
            inSites = true
        } else if (sequenceId && sessionId && 
                   line.startsWith(sessionId) && 
                   (inDomains || inSites)) {
            if (inDomains) {
                // #<session-ordinal>      <query-id[readingframe]>        <hit-type>      <PSSM-ID>       <from>  <to>    <E-Value>       <bitscore>      <accession>     <short-name>    <incomplete>    <superfamily PSSM-ID>
                def fields = line.split("\t")
                assert fields.size() == 12
                def hitType = fields[2]
                if (hitType.toUpperCase() == "SPECIFIC") {
                    def pssmId = fields[3]
                    def start = fields[4].toInteger()
                    def end = fields[5].toInteger()
                    def evalue = Double.parseDouble(fields[6])
                    def bitscore = Double.parseDouble(fields[7])
                    def modelAccession = fields[8]

                    // We need to use the PSSM ID to link domains and sites
                    def match
                    if (pssmHits.containsKey(pssmId)) {
                        match = pssmHits[pssmId]
                    } else {
                        def signature = new uk.ac.ebi.interpro.Signature(modelAccession, library)
                        match = new uk.ac.ebi.interpro.Match(modelAccession, signature)
                        pssmHits[pssmId] = match
                    }
                    match.addLocation(new uk.ac.ebi.interpro.Location(start, end, evalue, bitscore))
                }
            } else {
                // #<session-ordinal>      <query-id[readingframe]>        <annot-type>    <title> <residue(coordinates)>  <complete-size> <mapped-size>   <source-domain>
                def fields = line.split("\t")
                assert fields.size() == 8
                def hitType = fields[2]
                if (hitType.toUpperCase() == "SPECIFIC") {
                    def pssmId = fields[7]
                    if (pssmHits.containsKey(pssmId)) {
                        def description = fields[3]
                        def residues = pssmSites
                            .computeIfAbsent(pssmId) { [:] }
                            .computeIfAbsent(description) { [] as Set }
                        fields[4].split(",").each { residue ->
                            residues.add(residue)
                        }
                    }
                }
            }
        } else if (line.startsWith("ENDDOMAINS")) {
            assert inDomains == true
            assert inSites == false
            inDomains = false
        } else if (line.startsWith("ENDSITES")) {
            assert inDomains == false
            assert inSites == true
            pssmSites.each { pssmId, descriptionToResidues ->
                def match = pssmHits[pssmId]
                if (match != null) {
                    descriptionToResidues.each { description, residues ->
                        match.addSite(new uk.ac.ebi.interpro.Site(description, residues.join(",")))  // codenarc-disable-line JoinMismatchRule, JoinDuplicateRule
                    }
                }
            }
            pssmSites.clear()
            inSites = false
        } else if (line.startsWith("ENDQUERY")) {
            assert sequenceId != null
            
            def cddHits = [:]
            pssmHits.each { _key, match ->
                def modelAccession = match.modelAccession
                assert !cddHits.containsKey(modelAccession)
                cddHits[modelAccession] = match
            }
            
            pssmHits.clear()
            hits[sequenceId] = cddHits
            sequenceId = null
        } else if (line.startsWith("ENDSESSION")) {
            assert sessionId != null
            assert sequenceId == null
        } 
    }

    def filepath = task.workDir.resolve("cdd.json")
    filepath.text = groovy.json.JsonOutput.toJson(hits)
}
