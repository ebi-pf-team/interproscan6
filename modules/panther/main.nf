import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import java.nio.file.Files
import java.nio.file.Path
import uk.ac.ebi.interpro.FastaFile
import uk.ac.ebi.interpro.HMMER3
import uk.ac.ebi.interpro.Location
import uk.ac.ebi.interpro.Match
import uk.ac.ebi.interpro.TreeGrafter

process SEARCH_PANTHER {
    label 'mem_low', 'time_short', 'dynamic', 'ips6_container'

    input:
    tuple val(meta), val(meta2), path(fasta)
    path hmmdir
    val hmmfile

    output:
    tuple val(meta), val(meta2), path("hmmsearch.out")

    script:
    """
    hmmsearch \
        -Z 65000000 -E 0.001 \
        --cpu ${task.cpus} \
        ${hmmdir}/${hmmfile} ${fasta} > hmmsearch.out
    """
}

process PREPARE_TREEGRAFTER {
    label    'mem_low', 'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(meta2), val(hmmseach_out)
    val panther_dir
    val msf

    output:
    tuple val(meta), val(meta2), path("panther.json"),                             emit: json
    tuple val(meta), val(meta2), val(sequenceIds), val(familyIds), val(fastas),    emit: fasta
    
    exec:
    def hmmer_matches = HMMER3.parseOutput(hmmseach_out, "PANTHER")

    hmmer_matches = hmmer_matches.collectEntries { seqId, matches ->
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
                def familyId = (m1.modelAccession =~ /^(PTHR\d+)/)[0][1]
                m1.signature.accession = familyId
                def m2 = new Match(familyId, m1.evalue, m1.score, m1.bias, m1.signature)
                m2.included = m1.included
                // Only keep the domain with the highest score
                m2.locations = [locations.max { it.score }]
                // Init empty TreeGrafter attribute
                m2.treegrafter = new TreeGrafter(null)
                return m2
            }
            .findAll { it != null }

        def bestMatch = filteredMatches.max { it.score }
        return bestMatch ? [(seqId): [(bestMatch.modelAccession): bestMatch]] : [:]
    }

    def filepath = task.workDir.resolve("panther.json")
    filepath.text = JsonOutput.toJson(hmmer_matches)

    familyIds = []
    fastas = []
    sequenceIds = []
    hmmer_matches.each { seqId, matches ->
        // Ensure we only have one family
        assert matches.size() == 1
        def match = matches.values().first()
        
        // Ensure we only have one domain
        assert match.locations.size() == 1
        def location = match.locations.first()
        assert location.queryAlignment.length() == location.targetAlignment.length()

        // Get expected length of the sequence
        def familyId = match.modelAccession
        def fastaPath = file("${panther_dir}/${msf}/${familyId}.AN.fasta")
        assert fastaPath.exists()
        def sequences = FastaFile.parse(fastaPath)  // [md5 : "seq"]
        def length = sequences.values().first().length()

        // Query sequence to graft
        def sb = new StringBuilder()

        // Pad N-terminal
        sb << ("-" * (location.hmmStart - 1))

        // Build sequence
        def targetAlignment = location.targetAlignment.replaceAll(/(?i)[UO]/, 'X')
        for (int i = 0; i < targetAlignment.length(); i++) {
            def hmmChar = location.queryAlignment[i]
            def seqChar = targetAlignment[i]

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
        def sequence = sb.toString()
        assert sequence.length() == length

        def fastaFile = Files.createTempFile(task.workDir, "tmp", ".fa")
        fastaFile.withWriter { writer ->
            writer.writeLine(">${familyId}")
            for (int i = 0; i < sequence.length(); i += 80) {
                def j = Math.min(i + 80, length) - 1
                def line = sequence[i..j]
                writer.writeLine(line)                
            }
        }

        familyIds.add( familyId )
        fastas.add( fastaFile )
        sequenceIds.add( seqId )
    } 
}


process RUN_TREEGRAFTER {
    label 'mem_medium', 'time_short', 'dynamic', 'ips6_container'
    
    input:
    tuple val(meta), val(meta2), val(sequenceIds), val(familyIds), val(fastas)
    path panther_dir
    val msf

    output:
    tuple val(meta), val(meta2), path("epang.tsv")

    script:
    def commands = ""

    [sequenceIds, familyIds, fastas]
        .transpose()
        .each { entry -> 
            def seqID  = entry[0]
            def family = entry[1]
            def fastaPath = entry[2]
           
           // Run EPA-ng
            def epang_command = "epa-ng"
            epang_command += " -G 0.05"
            epang_command += " -m WAG"
            epang_command += " -T ${task.cpus}"
            epang_command += " -t ${panther_dir}/${msf}/${family}.bifurcate.newick"
            epang_command += " -s ${panther_dir}/${msf}/${family}.AN.fasta"
            epang_command += " -q ${fastaPath}"
            epang_command += " --redo"

            // Parse results
            def py_command = "python ${projectDir}/bin/panther/parse_epang.py"
            py_command += " epa_result.jplace"
            py_command += " ${panther_dir}/${msf}/${family}.newick"

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
    label    'mem_low', 'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(meta2), val(hmmseach_json), val(epagn_tsv)

    output:
    tuple val(meta), val(meta2), path("panther.json")

    exec:
    def matches = new JsonSlurper().parse(hmmseach_json).collectEntries { seqId, jsonMatches ->
        [(seqId): jsonMatches.collectEntries { matchId, jsonMatch ->
            [(matchId): Match.fromMap(jsonMatch)]
        }]
    }

    epagn_tsv.eachLine { line ->
        line = line.trim()
        def fields = line.split()
        assert fields.size() == 3
        def seqId = fields[0]
        def familyId = fields[1]
        def nodeId = fields[2]
        
        def match = matches[seqId]?.get(familyId)
        assert match != null
        match.treegrafter = new TreeGrafter(nodeId)
    }

    def filepath = task.workDir.resolve("panther.json")
    filepath.text = JsonOutput.toJson(matches)
}