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
    def matches = HMMER3.parseOutput(hmmsearch_out.toString())

    def matchesInfo = [:]
    JsonSlurper jsonSlurper = new JsonSlurper()
    def rules = jsonSlurper.parse(new File(rulesPath.toString()))
    matches.each { sequenceId, domains ->
        domains.each { modelId, domain ->
            def domHits = []
            def sortedLocations = domain.locations.sort { a, b ->
                a.evalue <=> b.evalue ?: b.score <=> a.score
            }
            sortedLocations.each { location ->
                def hmmFrom = location.hmmStart
                def hmmTo = location.hmmEnd
                def hmmAlign = location.querySequence
                def seqFrom = location.start
                def seqTo = location.end
                def seqAlign = location.targetSequence

                def map = mapHMMToSeq(hmmFrom, hmmAlign, seqAlign)

                def ruleSites = []
                def rule = rules.get(modelId, null)
                if (rule) {
                    rule.Groups.each { grp, positions ->
                        int passCount = 0
                        def positionsParsed = []
                        positions.eachWithIndex { pos, posNum ->
                            def condition = pos.condition.replaceAll(/[-()]/) { m ->
                                switch (m[0]) {
                                    case '-': ''; case '(': '{'; case ')': '}'
                                }
                            }.replace('x', '.')
                            println "positions: $positions"

                            def querySeq = seqAlign.replaceAll('-', '')
                            println "querySeq: $querySeq"
                            def targetSeq = (pos.hmmStart < map.size() && pos.hmmEnd < map.size()) ?
                                            querySeq[map[pos.hmmStart]..map[pos.hmmEnd]] : ''

                            if (targetSeq ==~ condition) {
                                passCount++
                                if (pos.start == 'Nter') pos.start = seqFrom
                                if (pos.end == 'Cter') pos.end = seqTo
                            }

                            positionsParsed << [
                                description: pos.desc,
                                group: pos.group as int,
                                hmmEnd: pos.hmmEnd,
                                hmmStart: pos.hmmStart,
                                label: pos.label,
                                numLocations: 1,
                                siteLocations: [[
                                    end: pos.end,
                                    residue: pos.condition,
                                    start: pos.start
                                ]]
                            ]
                        }
                        if (passCount == positions.size()) {
                            ruleSites.addAll(positionsParsed)
                        }
                    }
                }

                if (!ruleSites.isEmpty()) {
                    domHits << [
                        score: location.score as float,
                        evalue: location.evalue as float,
                        hmmStart: hmmFrom,
                        hmmEnd: hmmTo,
                        hmmAlign: hmmAlign,
                        start: seqFrom,
                        end: seqTo,
                        alignment: seqAlign,
                        sites: ruleSites,
                        hmmLength: location.hmmLength,
                        envelopeStart: location.envelopeStart as int,
                        envelopeEnd: location.envelopeEnd as int
                    ]
                }
            }

            domain.locations = domHits
            if (sortedLocations) {
                domain.score = sortedLocations[0].score as float
                domain.evalue = sortedLocations[0].evalue as float
            }
        }
    }

    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}


def mapHMMToSeq(int hmmPos, String hmm, String seq) {
    int seqPos = 0
    def map = [-1] * hmmPos + [0]

    hmm.eachWithIndex { character, i ->
        map[hmmPos++] = seqPos

        if (character != '.') hmmPos++
        if (seq[i] != '-') seqPos++
    }

    return map
}



