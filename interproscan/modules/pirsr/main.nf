import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import java.util.regex.Pattern

process RUN_PIRSR {
    label 'small'

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
    label 'small'

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

    def validMatches = [:]
    hmmerMatches.each { seqId, matches ->
        def filteredSeqMatches = [:]
        matches.each { modelAccession, match ->
            List<Location> sortedLocations = match.locations.sort { loc ->
                [loc.evalue, -loc.score]  // sorting by evalue ASC, score DESC
            }
            List<Location> selectedLocations = []
            sortedLocations.each { location ->
                if (!location.targetAlignment || !location.queryAlignment) {
                    return
                }
                List<Integer> map = mapHMMToSeq(location.hmmStart,
                                                location.queryAlignment,
                                                location.targetAlignment)

                List<Site> ruleSites = []
                def rule = rules.get(modelAccession, null)
                if (rule) {
                    rule.Groups.each { grp, positions ->
                        int passCount = 0
                        List<Site> positionsParsed = []
                        positions.each { pos ->
                            String condition = pos.condition.replaceAll(/[-()]/) { m ->
                                switch (m[0]) {
                                    case '-': return ''
                                    case '(': return '{'
                                    case ')': return '}'
                                }
                            }.replace('x', '.')

                            String querySeq = location.targetAlignment.replaceAll('-', '')
                            if (pos.hmmStart < map.size() && pos.hmmEnd < map.size()) {
                                targetSeq = querySeq[map[pos.hmmStart]..<map[pos.hmmEnd] + 1]
                                residue = location.targetAlignment[map[pos.hmmStart]..<map[pos.hmmEnd] + 1]
                            } else {
                                targetSeq = ''
                            }

                            if (targetSeq ==~ condition) {
                                passCount++
                                if (pos.start == 'Nter') pos.start = location.start
                                if (pos.end == 'Cter') pos.end = location.end
                            }
                            def (residueStart, residueEnd, residue) = [0, 0, null]
                            Map<Integer, Integer> seqAlignmentPosMap = getPositionMap(location.targetAlignment, location.start)
                            Map<Integer, Integer> seqAlignmentReversePosMap = seqAlignmentPosMap.collectEntries { k, v -> [(v): k] }
                            Map<Integer, Integer> hmmAlignmentPosMap = getPositionMap(location.queryAlignment, location.hmmStart)
                            if (hmmAlignmentPosMap.containsKey(pos.hmmStart)) {
                                int residueStartSeqAlign = hmmAlignmentPosMap[pos.hmmStart]
                                if (hmmAlignmentPosMap.containsKey(pos.hmmEnd)) {
                                    int residueEndSeqAlign = hmmAlignmentPosMap[pos.hmmEnd]
                                    residue = location.targetAlignment.substring(residueStartSeqAlign, residueEndSeqAlign + 1)
                                    if (seqAlignmentReversePosMap.containsKey(residueStartSeqAlign) &&
                                            seqAlignmentReversePosMap.containsKey(residueEndSeqAlign)) {
                                        residueStart = seqAlignmentReversePosMap[residueStartSeqAlign]
                                        residueEnd = seqAlignmentReversePosMap[residueEndSeqAlign]
                                    }
                                }
                            }
                            if (residueStart != 0 && residueEnd != 0) {
                                SiteLocation siteLocation = new SiteLocation(residue, residueStart, residueEnd)
                                positionsParsed << new Site(
                                    pos.desc,
                                    pos.group as int,
                                    pos.hmmEnd,
                                    pos.hmmStart,
                                    pos.label,
                                    [siteLocation]
                                )
                            }
                        }

                        if (passCount == positions.size()) {
                            ruleSites.addAll(positionsParsed)
                        }
                    }
                }
                if (!ruleSites.isEmpty()) {
                    location.sites = ruleSites
                    selectedLocations << location
                }
            }
            if (!selectedLocations.isEmpty()) {
                match.score = selectedLocations[0].score
                match.evalue = selectedLocations[0].evalue
                match.locations = selectedLocations
                filteredSeqMatches[modelAccession] = match
            }
        }
        if (!filteredSeqMatches.isEmpty()) {
            validMatches[seqId] = filteredSeqMatches
        }
    }
    def json = JsonOutput.toJson(validMatches)
    new File(outputFilePath.toString()).write(json)
}

def mapHMMToSeq(int hmmStart, String querySeq, String targetSeq) {
    int seqPos = 0
    def map = [0]
    for (int i = 1; i <= hmmStart; i++) {
        map << -1
    }
    for (int i = 0; i < querySeq.length(); i++) {
        map = map[0..hmmStart - 1] + [seqPos]
        if (querySeq[i] != '.') {
            hmmStart += 1
        }
        if (targetSeq[i] != '-') {
            seqPos += 1
        }
    }
    return map
}

def getPositionMap(alignment, position) {
    Map<Integer, Integer> positionMap = [:]
    alignment.eachWithIndex { character, index ->
        if (Character.isLetter(character as char)) {
            positionMap[position] = index
            position++
        }
    }
    return positionMap
}
