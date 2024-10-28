import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import java.util.regex.Pattern

process RUN_PIRSR {
    label 'hmmer_runner'

    input:
    tuple val(meta), path(fasta)
    path hmmdb

    output:
    tuple val(meta), path("hmmsearch.out")

    script:
    """
    /opt/hmmer3/bin/hmmsearch \
        -E 0.01 --acc \
        --cpu ${task.cpus} \
        ${hmmdb} ${fasta} > hmmsearch.out
    """
}

process PARSE_PIRSR {
    label 'analysis_parser'

    input:
    tuple val(meta), val(hmmsearch_out)
    val rulesPath

    output:
    tuple val(meta), path("pirsr.json")

    exec:
    def outputFilePath = task.workDir.resolve("pirsr.json")
    def hmmerMatches = HMMER3.parseOutput(hmmsearch_out.toString())

    JsonSlurper jsonSlurper = new JsonSlurper()
    def rules = jsonSlurper.parse(new File(rulesPath.toString()))

    hmmerMatches = hmmerMatches.collectEntries { seqId, matches ->
        matches.each { modelId, match ->
            def sortedLocations = match.locations.sort { a, b ->
                a.evalue <=> b.evalue ?: b.score <=> a.score
            }
            sortedLocations.each { location ->
                def map = mapHMMToSeq(location.hmmStart,
                                  location.querySequence,
                                  location.targetSequence)
//                 println "location.hmmStart: ${location.hmmStart}"
//                 println "Map: ${map}"

                def ruleSites = []
                def rule = rules.get(modelId, null)
                if (rule) {
                    rule.Groups.each { grp, positions ->
                        int passCount = 0
                        def positionsParsed = []
                        positions.each { pos ->
                            def condition = pos.condition.replaceAll(/[-()]/) { m ->
                                switch (m[0]) {
                                    case '-': ''; case '(': '{'; case ')': '}'
                                }
                            }.replace('x', '.')
                            def querySeq = location.targetSequence.replaceAll('-', '')

                            if (pos.hmmStart < map.size() && pos.hmmEnd < map.size()) {
                                targetSeq = querySeq[map[pos.hmmStart]..map[pos.hmmEnd]]
                            } else {
                                targetSeq = ''
                            }

                            if (targetSeq ==~ condition) {
                                passCount++
                                if (pos.start == 'Nter') pos.start = location.start
                                if (pos.end == 'Cter') pos.end = location.end
                            }

                            SiteLocation siteLocation = new SiteLocation(
                                pos.condition, pos.start, pos.end)

                            positionsParsed << [new Site(
                                pos.desc,
                                pos.group as int,
                                pos.hmmEnd,
                                pos.hmmStart,
                                pos.label,
                                [siteLocation]
                            )]
                        }
                        if (passCount == positions.size()) {
                            ruleSites.addAll(positionsParsed)
                        }
                    }
                }

                if (!ruleSites.isEmpty()) {
                    location.sites = ruleSites
                }
            }

            if (sortedLocations) {
                match.score = sortedLocations[0].score
                match.evalue = sortedLocations[0].evalue
            }
        }
    }

    def json = JsonOutput.toJson(hmmerMatches)
    new File(outputFilePath.toString()).write(json)
}

def mapHMMToSeq(int hmmStart, String querySeq, String targetSeq) {
    // map base positions from alignment, from query HMM coords to (ungapped) target sequence coords
    int seqPos = 0
    def map = [0] + (1..<hmmStart).collect { -1 }
    querySeq.eachWithIndex { character, i ->
        map << seqPos
        if (character != '.') hmmStart++
        if (targetSeq[i] != '-') seqPos++
    }
    return map
}
