process PARSE_PIRSR {
    label    'mem_low', 'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(hmmsearch_out)
    val rulesfile

    output:
    tuple val(meta), path("pirsr.json")

    exec:
    def hmmerMatches = uk.ac.ebi.interpro.HMMER3.parseOutput(hmmsearch_out, "PIRSR")
    def rules = new groovy.json.JsonSlurper().parse(rulesfile)
    def validMatches = [:]
    hmmerMatches.each { seqId, matches ->
        def filteredSeqMatches = [:]
        matches.each { modelAccession, match ->
            // set the signature name, which defaults to null in the HMMER3 parser
            match.signature.name = modelAccession

            def sortedLocations = match.locations.sort { loc ->
                [loc.evalue, -loc.score]  // sorting by evalue ASC, score DESC
            }
            def selectedLocations = []
            sortedLocations.each { location ->
                if (!location.included || !location.targetAlignment || !location.queryAlignment) {
                    return
                }
                def map = mapHMMToSeq(location.hmmStart,
                                                location.queryAlignment,
                                                location.targetAlignment)

                def ruleSites = [] as Set
                def rule = rules.get(modelAccession, null)
                if (rule) {
                    rule.Groups.each { grp, positions ->
                        def passCount = 0
                        def positionsParsed = [] as Set
                        positions.each { pos ->
                            def condition = pos.condition.replaceAll(/[-()]/) { m ->
                                if (m[0] == '-') {
                                    return ''
                                } else if (m[0] == '(') {
                                    return '{'
                                } else if (m[0] == ')') {
                                    return '}'
                                }
                            }.replace('x', '.')

                            def querySeq = location.targetAlignment.replaceAll('-', '')
                            def targetSeq
                            if (pos.hmmStart < map.size() && pos.hmmEnd < map.size()) {
                                targetSeq = querySeq[map[pos.hmmStart]..<map[pos.hmmEnd] + 1]
                            } else {
                                targetSeq = ''
                            }

                            if (targetSeq ==~ condition) {
                                passCount += 1
                                if (pos.start == 'Nter') pos.start = location.start
                                if (pos.end == 'Cter') pos.end = location.end
                            }
                            def (residueStart, residueEnd, residue) = [0, 0, null]
                            def seqAlignmentPosMap = getPositionMap(location.targetAlignment, location.start)
                            def seqAlignmentReversePosMap = seqAlignmentPosMap.collectEntries { k, v -> [(v): k] }
                            def hmmAlignmentPosMap = getPositionMap(location.queryAlignment, location.hmmStart)
                            if (hmmAlignmentPosMap.containsKey(pos.hmmStart)) {
                                def residueStartSeqAlign = hmmAlignmentPosMap[pos.hmmStart]
                                if (hmmAlignmentPosMap.containsKey(pos.hmmEnd)) {
                                    def residueEndSeqAlign = hmmAlignmentPosMap[pos.hmmEnd]
                                    residue = location.targetAlignment.substring(residueStartSeqAlign, residueEndSeqAlign + 1)
                                    if (seqAlignmentReversePosMap.containsKey(residueStartSeqAlign) &&
                                            seqAlignmentReversePosMap.containsKey(residueEndSeqAlign)) {
                                        residueStart = seqAlignmentReversePosMap[residueStartSeqAlign]
                                        residueEnd = seqAlignmentReversePosMap[residueEndSeqAlign]
                                    }
                                }
                            }
                            if (residueStart != 0 && residueEnd != 0) {
                                positionsParsed << [
                                    pos.desc,
                                    residue,
                                    residueStart,
                                    residueEnd
                                ]
                            }
                        }

                        if (passCount == positions.size()) {
                            ruleSites.addAll(positionsParsed)
                        }
                    }
                }
                if (!ruleSites.isEmpty()) {
                    def groupedPositions = [:]
                    ruleSites.each { desc, residue, residueStart, residueEnd ->
                        def siteLocations = groupedPositions.computeIfAbsent(desc) { [] }
                        siteLocations << new uk.ac.ebi.interpro.SiteLocation(residue, residueStart, residueEnd)
                    }
                    groupedPositions.each { description, siteLocations ->
                        siteLocations.sort { a, b -> a.start <=> b.start ?: a.end <=> b.end }
                        location.addSite(new uk.ac.ebi.interpro.Site(description, siteLocations))
                    }
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

    def filepath = task.workDir.resolve("pirsr.json")
    filepath.text  = groovy.json.JsonOutput.toJson(validMatches)
}

def mapHMMToSeq(int hmmStart, String querySeq, String targetSeq) {
    def seqPos = 0
    def currentHmmStart = hmmStart
    def map = ([0] + (1..hmmStart).collect { -1 }) as List<Integer>

    querySeq.eachWithIndex { char queryChar, int index ->
        map[currentHmmStart] = seqPos
        if (queryChar != '.') {
            currentHmmStart += 1
            map << -1
        }
        if (targetSeq[index] != '-') {
            seqPos += 1
        }
    }

    return map
}

def getPositionMap(alignment, position) {
    def positionMap = [:]
    alignment.eachWithIndex { character, index ->
        if (Character.isLetter(character as char)) {
            positionMap[position] = index
            position += 1
        }
    }
    return positionMap
}
