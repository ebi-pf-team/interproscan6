import com.fasterxml.jackson.databind.ObjectMapper
import groovy.json.JsonOutput
import uk.ac.ebi.interpro.FastaFile

process REPORT_NO_MATCHES {
    label    'mem_low', 'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(fasta), val(member_matches)

    output:
    tuple val(meta), path("no_matches.json")

    exec:
    def seqsWithMatches = [] as Set
    member_matches.each { matchesPath ->
        def matchesFileMap = matchesPath.newReader().withCloseable { reader ->
            new ObjectMapper().readValue(reader, Map)
        }
        matchesFileMap.each { String seqMd5, Map matches ->
            seqsWithMatches.add(seqMd5)
        }
    }

    def sequences = FastaFile.parse(fasta)  // [md5: sequence]
    def noMatches = [:]
    sequences.each { String seqMd5, String seq ->
        if (!seqsWithMatches.contains(seqMd5)) {
            noMatches[seqMd5] = [:]
        }
    }
    def output_file = task.workDir.resolve("no_matches.json")
    output_file.text = JsonOutput.toJson(noMatches)
}
