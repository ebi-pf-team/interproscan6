import groovy.json.JsonOutput
import groovy.json.JsonSlurper

process PFSCAN_RUNNER {
    /*
    The ps_scan.pl script is a wrapper for the pfscan tool that is provided by the
    pftools developers. It automates running pfscan for all provided patterns and
    includes post-processing of the hits.
    */
    label 'prosite_pfscan_runner'

    input:
        tuple val(meta), path(fasta)
        path data
        path evaluator

    output:
        tuple val(meta), path("ps_scan.out")

    script:
    """
        perl /opt/pftools/ps_scan.pl \
        ${fasta} \
        -d ${data} \
        --pfscan /opt/pftools/pfscanV3 \
        -b ${evaluator} \
        -r -s -o ipro > ps_scan.out
    """
}


process PFSCAN_PARSER {
    label 'analysis_parser'

    input:
        tuple val(meta), val(pfscan_out)

    output:
        tuple val(meta), path("ps_scan_parsed.json")

    exec:
    Map patternsMatches = [:]
    pfscan_out.eachLine { line ->
        String lineStrip = line.trim()
        if (!lineStrip || lineStrip.startsWith("pfscanV3 is not meant to be used with a single profile")) {
            return
        }
        List<String> matchInfo = lineStrip.split('\t')
        if (matchInfo.size() < 9) {
            return
        }
        List<String> matchDetails = matchInfo[8].split(';')
        String level = matchDetails[1].trim()
        if (!level.startsWith("LevelTag") || !level.contains("0")) {  // apenas matches fortes
            return
        }

        String seqId = matchInfo[0]
        String matchId = matchInfo[2]
        String alignment = matchDetails[2].replaceAll('Sequence ', '').replaceAll('"', '').replaceAll('\\.', '').trim()
//         List cigarAlignment = CigarAlignmentParser.parse(alignment)

        start = matchInfo[3].toInteger()
        end = matchInfo[4].toInteger()
        level = "STRONG"
        alignment = alignment
        cigarAlignment = "" // CigarAlignmentParser.encode(cigarAlignment),

        if (patternsMatches.containsKey(seqId)) {
            match = patternsMatches[seqId]
        } else {
            match = new Match(matchId)
            patternsMatches[seqId] = match
        }
        Location location = new Location(start, end, level, alignment, cigarAlignment)
        match.addLocation(location)
    }

    def outputFilePath = task.workDir.resolve("ps_scan_parsed.json")
    def json = JsonOutput.toJson(patternsMatches)
    new File(outputFilePath.toString()).write(json)
}
