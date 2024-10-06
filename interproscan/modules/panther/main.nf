import groovy.json.JsonOutput
import groovy.json.JsonSlurper

process SEARCH_PANTHER {
    label 'hmmer_runner'

    input:
    tuple val(meta), path(fasta)
    path hmmdb

    output:
    tuple val(meta), path("hmmsearch.out")

    script:
    """
    /opt/hmmer3/bin/hmmsearch \
        -Z 65000000 -E 0.001 --domE 0.00000001 --incdomE 0.00000001 \
        --cpu ${task.cpus} \
        ${hmmdb} ${fasta} > hmmsearch.out
    """
}

process PREPARE_TREEGRAFTER {
    input:
    tuple val(meta), val(hmmseach_out)
    val msf_dir

    output:
    tuple val(meta), path("panther.json"),                             emit: json
    tuple val(meta), val(sequenceIds), val(familyIds), path(fastas),   emit: fasta
    
    exec:
    def hmmerMatches = HMMER3.parseOutput(hmmseach_out.toString())

    hmmerMatches = hmmerMatches.collectEntries { seqId, matches ->
        // Filter matches to only those with locations that have a score > 100
        def filteredMatches = matches
            .values()
            .collect { m1 ->
                // TODO: do we need to filter by e-value here? Cannot it be done with hmmsearch?
                if (!m1.included || m1.evalue > 0.00000001) {
                    return null
                }

                def locations = m1.locations.findAll { it.included }
                if (!locations) {
                    return null
                }

                // Rename model accession (PTHR23076.orig.30.pir -> PTHR23076)
                String familyId = (m1.modelAccession =~ /^(PTHR\d+)/)[0][1]

                Match m2 = new Match(familyId, m1.evalue, m1.score, m1.bias)
                m2.included = m1.included
                // Only keep the domain with the highest score
                m2.locations = [locations.max { it.score }]
                return m2
            }
            .findAll { it != null }

        Match bestMatch = filteredMatches.max { it.locations[0].score }
        return bestMatch ? [(seqId): [(bestMatch.modelAccession): bestMatch]] : [:]
    }

    def outputFilePath = task.workDir.resolve("panther.json")
    def json = JsonOutput.toJson(hmmerMatches)
    new File(outputFilePath.toString()).write(json)

    familyIds = []
    fastas = []
    sequenceIds = []
    String msfPath = msf_dir.toString()
    hmmerMatches.each { seqId, matches ->
        // Ensure we only have one family
        assert matches.size() == 1
        Match match = matches.values().first()
        
        // Ensure we only have one domain
        assert match.locations.size() == 1
        Location location = match.locations.first()
        assert location.querySequence.length() == location.targetSequence.length()

        // Get expected length of the sequence
        String familyId = match.modelAccession
        Path fastaPath = file("${msfPath}/${familyId}.AN.fasta")
        assert fastaPath.exists()
        def sequences = FastaFile.parse(fastaPath.toString())
        int length = sequences.first().sequence.length()

        // Query sequence to graft
        StringBuilder sb = new StringBuilder()

        // Pad N-terminal
        sb << ("-" * (location.hmmStart - 1))

        // Build sequence
        String targetSequence = location.targetSequence.replaceAll(/(?i)[UO]/, 'X')
        for (int i = 0; i < targetSequence.length(); i++) {
            char hmmChar = location.querySequence[i]
            char seqChar = targetSequence[i]

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

        fastaPath = task.workDir.resolve("${seqId}.fa")
        File file = new File(fastaPath.toString())
        file.withWriter { writer ->
            writer.writeLine(">${familyId}")
            for (int i = 0; i < sequence.length(); i += 80) {
                int j = Math.min(i + 80, length) - 1
                String line = sequence[i..j]
                writer.writeLine(line)                
            }
        }

        familyIds.add( familyId )
        fastas.add( fastaPath )
        sequenceIds.add( seqId )
    } 
}


process RUN_TREEGRAFTER {
    label 'analysis_parser'
    
    input:
    tuple val(meta), val(sequenceIds), val(familyIds), path(fastas)
    path msf_dir

    output:
    tuple val(meta), path("epang.tsv")

    script:
    def commands = ""

    [sequenceIds, familyIds, fastas]
        .transpose()
        .eachWithIndex { entry, index -> 
            String seqID  = entry[0]
            String family = entry[1]
            def fastaPath = entry[2]
           
           // Run EPA-ng
            def command = "/opt/epa-ng/bin/epa-ng"
            command += " -G 0.05"
            command += " -m WAG"
            command += " -T ${task.cpus}"
            command += " -t ${msf_dir.toString()}/${family}.bifurcate.newick"
            command += " -s ${msf_dir.toString()}/${family}.AN.fasta"
            command += " -q ${fastaPath}"
            command += " --redo\n"

            // Parse results
            command += "python ${projectDir}/bin/panther/parse_epang.py"
            command += " epa_result.jplace"
            command += " ${msf_dir.toString()}/${family}.newick | "

            // Add sequence ID
            command += "awk '{print \"${seqID}\",\$0}' >> epang.tsv\n"

            commands += command 
        }

    """
    ${commands}
    """
}

process PARSE_PANTHER {
    label 'analysis_parser'

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