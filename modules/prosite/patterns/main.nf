import groovy.json.JsonOutput
import groovy.json.JsonSlurper

process RUN_PFSCAN {
    /*
    The ps_scan.pl script is a wrapper for the pfscan tool that is provided by the
    pftools developers. It automates running pfscan for all provided patterns and
    includes post-processing of the hits.
    */
    label 'mini', 'ips6_container'

    input:
        tuple val(meta), path(fasta)
        path dirpath
        val datfile
        val evafile

    output:
        tuple val(meta), path("ps_scan.out")

    script:
    """
        perl ${projectDir}/bin/prosite/ps_scan.pl \
        ${fasta} \
        -d ${dirpath}/${datfile} \
        --pfscan pfscanV3 \
        -b ${dirpath}/${evafile} \
        -r -s -o ipro > ps_scan.out
    """
}


process PARSE_PFSCAN {
    label    'tiny'
    executor 'local'

    input:
        tuple val(meta), val(pfscan_out)

    output:
        tuple val(meta), path("prositepatterns.json")

    exec:
    Map<String, Map<String, Match>> patternsMatches = [:]
    SignatureLibraryRelease library = new SignatureLibraryRelease("PROSITE patterns", null)
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
            return // skipping non-strong matches
        } else {
            level = "STRONG"
        }
        String alignment = matchDetails[2].replaceAll('Sequence ', '').replaceAll('"', '').replaceAll('\\.', '').trim()
        String cigarAlignment = Match.encodeCigarAlignment(alignment)
        patternsMatches.computeIfAbsent(seqId) { [:] }
        Match matchObj = patternsMatches[seqId].computeIfAbsent(modelAccession) {
            new Match(modelAccession, new Signature(modelAccession, library))
        }
        Location location = new Location(start, end, level, alignment, cigarAlignment)
        matchObj.addLocation(location)
    }

    def outputFilePath = task.workDir.resolve("prositepatterns.json")
    def json = JsonOutput.toJson(patternsMatches)
    new File(outputFilePath.toString()).write(json)
}
