import groovy.json.JsonOutput

process RUN_RPSBLAST {
    label 'small', 'ips6_container'

    input:
    tuple val(meta), path(fasta)
    val library

    output:
    tuple val(meta), path("rpsblast.out")

    script:
    """
    export LD_LIBRARY_PATH="/opt/blast/lib"
    /opt/blast/rpsblast \
        -query ${fasta} \
        -db "${library}/Cdd_NCBI" \
        -out rpsblast.out \
        -evalue 0.01 -seg no -outfmt 11
    """
}


process RUN_RPSPROC {
    /*
    Process output from offline rpsblast to annotate sequence with conserved
    domain information

    rpsbproc desc:
    This utility processes domain hit data produced by local RPS-BLAST and
    generates domain family and/or superfamily annotations on the query
    sequences. Instead of retrieving domain data from database, this program
    processes dumped datafiles to obtain required information. All data files
    are downloadable from NCBI ftp site. Read README file for details
    */
    label 'small', 'ips6_container'

    input:
    tuple val(meta), val(rpsblast_out)
    path library

    output:
    tuple val(meta), path("rpsbproc.out")

    script:
    """
    /opt/rpsbproc/RpsbProc-x64-linux/rpsbproc \
        --infile ${rpsblast_out} \
        --outfile rpsbproc.out \
        --data-path ${library} \
        -m std
    """
}


process PARSE_RPSPROC {
    label 'small'

    input:
    tuple val(meta), val(rpsbproc_out)

    output:
    tuple val(meta), path("cdd.json")

    exec:
    String sessionId = null
    String sequenceId = null
    boolean inDomains = false
    boolean inSites = false
    def hits = [:]
    def pssmHits = [:]
    file(rpsbproc_out.toString()).eachLine { line -> 
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
                String hitType = fields[2]
                if (hitType.toUpperCase() == "SPECIFIC") {
                    String pssmId = fields[3]
                    int start = fields[4].toInteger()
                    int end = fields[5].toInteger()
                    Double evalue = Double.parseDouble(fields[6])
                    Double bitscore = Double.parseDouble(fields[7])
                    String modelAccession = fields[8]

                    // We need to use the PSSM ID to link domains and sites
                    Match match
                    if (pssmHits.containsKey(pssmId)) {
                        match = pssmHits[pssmId]
                    } else {
                        match = new Match(modelAccession, evalue, bitscore)
                        pssmHits[pssmId] = match
                    }

                    match.addLocation(new Location(start, end))
                }
            } else {
                // #<session-ordinal>      <query-id[readingframe]>        <annot-type>    <title> <residue(coordinates)>  <complete-size> <mapped-size>   <source-domain>
                def fields = line.split("\t")
                assert fields.size() == 8
                String hitType = fields[2]
                if (hitType.toUpperCase() == "SPECIFIC") {
                    String pssmId = fields[7]
                    if (pssmHits.containsKey(pssmId)) {
                        String description = fields[3]
                        String residues = fields[4]
                        Site site = new Site(description, residues)
                        Match match = pssmHits[pssmId]
                        match.addSite(site)                        
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
            inSites = false
        } else if (line.startsWith("ENDQUERY")) {
            assert sequenceId != null
            
            def cddHits = [:]
            pssmHits.each { key, match ->
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

    def outputFilePath = task.workDir.resolve("cdd.json")
    def json = JsonOutput.toJson(hits)
    new File(outputFilePath.toString()).write(json)
}
