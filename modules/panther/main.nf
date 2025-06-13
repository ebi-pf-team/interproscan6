import groovy.json.JsonOutput
import groovy.json.JsonSlurper

process PREPARE_TREEGRAFTER {
    label    'tiny'
    executor 'local'

    input:
    tuple val(meta), val(hmmseach_out)
    val dir
    val msf

    output:
    tuple val(meta), path("panther.json"),                             emit: json
    tuple val(meta), val(sequenceIds), val(familyIds), val(fastas),    emit: fasta
    
    exec:
    def hmmerMatches = HMMER3.parseOutput(hmmseach_out.toString(), "PANTHER")

    hmmerMatches = hmmerMatches.collectEntries { seqId, matches ->
        // Filter matches to only those with locations that have a score > 100
        def filteredMatches = matches
            .values()
            .collect { m1 ->
                if (!m1.included) {
                    return null
                }

                def locations = m1.locations.findAll { it.included }
                if (!locations) {
                    return null
                }

                // Rename model accession (PTHR23076.orig.30.pir -> PTHR23076)
                String familyId = (m1.modelAccession =~ /^(PTHR\d+)/)[0][1]
                m1.signature.accession = familyId
                Match m2 = new Match(familyId, m1.evalue, m1.score, m1.bias, m1.signature)
                m2.included = m1.included
                // Only keep the domain with the highest score
                m2.locations = [locations.max { it.score }]
                // Init empty TreeGrafter attribute
                m2.treegrafter = new TreeGrafter(null)
                return m2
            }
            .findAll { it != null }

        Match bestMatch = filteredMatches.max { it.score }
        return bestMatch ? [(seqId): [(bestMatch.modelAccession): bestMatch]] : [:]
    }

    def outputFilePath = task.workDir.resolve("panther.json")
    def json = JsonOutput.toJson(hmmerMatches)
    new File(outputFilePath.toString()).write(json)

    familyIds = []
    fastas = []
    sequenceIds = []
    String dirPath = dir.toString()
    hmmerMatches.each { seqId, matches ->
        // Ensure we only have one family
        assert matches.size() == 1
        Match match = matches.values().first()
        
        // Ensure we only have one domain
        assert match.locations.size() == 1
        Location location = match.locations.first()
        assert location.queryAlignment.length() == location.targetAlignment.length()

        // Get expected length of the sequence
        String familyId = match.modelAccession
        Path fastaPath = file("${dirPath}/${msf}/${familyId}.AN.fasta")
        assert fastaPath.exists()
        Map sequences = FastaFile.parse(fastaPath.toString())  // [md5 : "seq"]
        int length = sequences.values().first().length()

        // Query sequence to graft
        StringBuilder sb = new StringBuilder()

        // Pad N-terminal
        sb << ("-" * (location.hmmStart - 1))

        // Build sequence
        String targetAlignment = location.targetAlignment.replaceAll(/(?i)[UO]/, 'X')
        for (int i = 0; i < targetAlignment.length(); i++) {
            char hmmChar = location.queryAlignment[i]
            char seqChar = targetAlignment[i]

            if (hmmChar != '.') {
                sb << seqChar
            }
        }

        // Pad C-terminal
        assert sb.length() <= length
        while (sb.length() < length) {
            sb << "-"
        }

        assert sb.length() <= length
        String sequence = sb.toString()
        assert sequence.length() == length

        File fastaFile = File.createTempFile("tmp", ".fa", task.workDir.toFile())
        fastaFile.withWriter { writer ->
            writer.writeLine(">${familyId}")
            for (int i = 0; i < sequence.length(); i += 80) {
                int j = Math.min(i + 80, length) - 1
                String line = sequence[i..j]
                writer.writeLine(line)                
            }
        }

        familyIds.add( familyId )
        fastas.add( fastaFile.toString() )
        sequenceIds.add( seqId )
    } 
}


process RUN_TREEGRAFTER {
    label 'small', 'ips6_container'
    
    input:
    tuple val(meta), val(sequenceIds), val(familyIds), val(fastas)
    path dir
    val msf

    output:
    tuple val(meta), path("epang.tsv")

    script:
    def commands = ""

    [sequenceIds, familyIds, fastas]
        .transpose()
        .each { entry -> 
            String seqID  = entry[0]
            String family = entry[1]
            def fastaPath = entry[2]
           
           // Run EPA-ng
            def epang_command = "/opt/epa-ng/bin/epa-ng"
            epang_command += " -G 0.05"
            epang_command += " -m WAG"
            epang_command += " -T ${task.cpus}"
            epang_command += " -t ${dir.toString()}/${msf}/${family}.bifurcate.newick"
            epang_command += " -s ${dir.toString()}/${msf}/${family}.AN.fasta"
            epang_command += " -q ${fastaPath}"
            epang_command += " --redo"

            // Parse results
            def py_command = "python ${projectDir}/bin/panther/parse_epang.py"
            py_command += " epa_result.jplace"
            py_command += " ${dir.toString()}/${msf}/${family}.newick"

            // Add sequence ID
            def awk_command = "awk '{print \"${seqID}\",\$0}'"

            // Only run Python + Awk if EPA-ng doesn't fail
            commands += "{ ${epang_command} && ${py_command} | ${awk_command} >> epang.tsv; } || :\n"
        }

    """
    touch epang.tsv
    ${commands}
    """
}

process PARSE_PANTHER {
    label    'tiny'
    executor 'local'

    input:
    tuple val(meta), val(hmmseach_json), val(epagn_tsv)

    output:
    tuple val(meta), path("panther.json")

    exec:
    File jsonFile = new File(hmmseach_json.toString())
    JsonSlurper jsonSlurper = new JsonSlurper()
    def matches = jsonSlurper.parse(jsonFile).collectEntries { seqId, jsonMatches ->
        [(seqId): jsonMatches.collectEntries { matchId, jsonMatch ->
            [(matchId): Match.fromMap(jsonMatch)]
        }]
    }

    file(epagn_tsv.toString()).eachLine { line ->
        line = line.trim()
        def fields = line.split()
        assert fields.size() == 3
        String seqId = fields[0]
        String familyId = fields[1]
        String nodeId = fields[2]
        
        Match match = matches[seqId]?.get(familyId)
        assert match != null
        match.treegrafter = new TreeGrafter(nodeId)
    }

    def outputFilePath = task.workDir.resolve("panther.json")
    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}