import groovy.json.JsonOutput
import groovy.json.JsonSlurper

process RUN_PFSCAN {
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


process PARSE_PFSCAN {
    label 'analysis_parser'

    input:
        tuple val(meta), val(pfscan_out)

    output:
        tuple val(meta), path("psscan.json")

    exec:
    Map<String, Map<String, Match>> patternsMatches = [:]
    pfscan_out.eachLine { line ->
        line = line.trim()
        if (!line || line.startsWith("pfscanV3 is not meant to be used with a single profile")) {
            return
        }
        List<String> matchInfo = line.split('\t')
        if (matchInfo.size() < 9) {
            return
        }
        String seqId = matchInfo[0]
        String modelAccession = matchInfo[2]
        int start = matchInfo[3].toInteger()
        int end = matchInfo[4].toInteger()

        List<String> matchDetails = matchInfo[8].split(';')
        String name = matchDetails[0].trim()
        String level = matchDetails[1].trim()
        if (!level.startsWith("LevelTag") || !level.contains("0")) {
            // skipping not strong matches
        } else {
            level = "STRONG"
        }
        String alignment = matchDetails[2].replaceAll('Sequence ', '').replaceAll('"', '').replaceAll('\\.', '').trim()
        String cigarAlignment = parseCigarAlignment(alignment)
        cigarAlignment = encodeCigarAlignment(cigarAlignment)
        if (patternsMatches.containsKey(seqId)) {
            match = patternsMatches[seqId]
        } else {
            match = new Match(modelAccession)
            patternsMatches[seqId] = match
        }
        Location location = new Location(start, end, level, alignment, cigarAlignment)
        match.addLocation(location)
    }

    def outputFilePath = task.workDir.resolve("psscan.json")
    def json = JsonOutput.toJson(patternsMatches)
    new File(outputFilePath.toString()).write(json)
}

def parseCigarAlignment(String alignment) {
    String cigarAlignment = ""
    alignment.each { baseChar ->
        if (baseChar.isUpperCase()) {
            cigarAlignment += "M" // match char
        } else if (baseChar.isLowerCase()) {
            cigarAlignment += "I" // insert char
        } else if (baseChar == "-") {
            cigarAlignment += "D" // delete char
        } else {
            throw new IllegalArgumentException("Unrecognized character ${baseChar} in ${alignment}")
        }
    }
    return cigarAlignment
}

def encodeCigarAlignment(String cigarAlignment) {  // Compress alignment, to give '5M' instead of 'MMMMM'
    if (!cigarAlignment) return ""
    cigarAlignment
        .split('')
        .inject([]) { groups, alignChar ->
            if (groups && groups.last()[1] == alignChar) {
                groups.last()[0]++  // cigar alignment char matches previous so add 1 to the count
            } else {
                groups << [1, alignChar]  // change type of alignment char, restart count
            }
            groups
        }
        .collect { count, alignChar -> "${count}${alignChar}" }
        .join()
}
