import com.fasterxml.jackson.databind.ObjectMapper
import groovy.json.JsonOutput

process REPORT_NO_MATCHES {
    executor 'local'

    input:
    tuple val(meta), val(member_matches)
    tuple val(meta), val(fasta)

    output:
    tuple val(meta), path("no_matches.json")

    exec:
    def seqsWithMatches = [] as Set
    member_matches.each { matchesPath ->
        def matchesFileMap = new ObjectMapper().readValue(new File(matchesPath.toString()), Map)
        matchesFileMap.each { String seqMd5, Map matches ->
            seqsWithMatches.add(seqMd5)
        }
    }

    Map<String, String> sequences = FastaFile.parse(fasta.toString())  // [md5: sequence]
    String outputFilePath = task.workDir.resolve("no_matches.json")
    def noMatches = [:]
    sequences.each { String seqMd5, String seq ->
        if (!seqsWithMatches.contains(seqMd5)) {
            noMatches[seqMd5] = [:]
        }
    }
    def json = JsonOutput.toJson(noMatches)
    new File(outputFilePath.toString()).write(json)
}
