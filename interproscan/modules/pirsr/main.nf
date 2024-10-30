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
        matches.each { modelAccession, match ->
            def sortedLocations = match.locations.sort { loc ->
                [loc.evalue, -loc.score]  // sorting by evalue ASC, score DESC
            }
            sortedLocations.each { location ->
                def map = mapHMMToSeq(location.hmmStart,
                                  location.querySequence,
                                  location.targetSequence)

                def ruleSites = []
                def rule = rules.get(modelAccession, null)
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

                            String targetSeq = ''
                            if (pos.hmmStart < map.size() && pos.hmmEnd < map.size()) {
                                targetSeq = querySeq[map[pos.hmmStart]..map[pos.hmmEnd]]
                            }

                            if (targetSeq ==~ condition) {
                                passCount++
                                if (pos.start == 'Nter') pos.start = location.start
                                if (pos.end == 'Cter') pos.end = location.end
                            }

                            def (residueStart, residueEnd, residue) = [0, 0, null]
                            def seqAlignmentPosMap = getPositionMap(location.targetSequence, location.start)
                            def seqAlignmentReversePosMap = seqAlignmentPosMap.collectEntries { k, v -> [(v): k] }
                            def hmmAlignmentPosMap = getPositionMap(location.querySequence, pos.hmmStart)
                            if (pos.hmmStart in hmmAlignmentPosMap) {
                                int residueStartSeqAlign = hmmAlignmentPosMap[pos.hmmStart]
                                if (pos.hmmEnd in hmmAlignmentPosMap) {
                                    int residueEndSeqAlign = hmmAlignmentPosMap[pos.hmmEnd]
                                    residue = location.targetSequence.substring(residueStartSeqAlign, residueEndSeqAlign + 1)
                                    if (residueStartSeqAlign in seqAlignmentReversePosMap &&
                                            residueEndSeqAlign in seqAlignmentReversePosMap) {
                                        residueStart = seqAlignmentReversePosMap[residueStartSeqAlign]
                                        residueEnd = seqAlignmentReversePosMap[residueEndSeqAlign]
                                    }
                                }
                            }

                            if (!residueStart == 0 && !residueEnd == 0) {
                                SiteLocation siteLocation = new SiteLocation(pos.condition, residueStart, residueEnd)
                                positionsParsed << [new Site(
                                    pos.desc,
                                    pos.group as int,
                                    pos.hmmEnd,
                                    pos.hmmStart,
                                    pos.label,
                                    [siteLocation]
                                )]
                            }
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
        return [(seqId): matches]
    }

    def json = JsonOutput.toJson(hmmerMatches)
    new File(outputFilePath.toString()).write(json)
}

def mapHMMToSeq(int hmmStart, String querySeq, String targetSeq) {
    int seqPos = 0
    def map = (0..<hmmStart).collect { -1 }
    querySeq.eachWithIndex { character, i ->
        map[hmmStart] = seqPos
        hmmStart += (character != '.') ? 1 : 0
        seqPos += (targetSeq[i] != '-') ? 1 : 0
    }
    return map
}

def getPositionMap(alignment, position) {
    def positionMap = [:]
    alignment.eachWithIndex { character, index ->
        if (Character.isLetter(character as char)) {
            positionMap[position] = index
            position++
        }
    }
    return positionMap
}
