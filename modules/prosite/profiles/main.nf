import groovy.io.FileType
import groovy.json.JsonOutput

import Match

process RUN_PFSEARCH {
    label 'tiny', 'ips6_container'

    input:
        tuple val(meta), path(fasta)
        path dirpath
        val profiles_dir

    output:
        tuple val(meta), path("prosite_profiles.out")

    script:
    """
    touch prosite_profiles.out
    find ${dirpath}/${profiles_dir} -type f | while read profile; do
        output=\$(/opt/pftools/pfsearchV3 "\${profile}" "${fasta}" -f -o 7 -t 4)
        if [[ -n "\${output}" ]]; then
            echo "\${output}" >> prosite_profiles.out
        fi
    done
    """
}

process PARSE_PFSEARCH {
    label    'tiny'
    executor 'local'

    input:
        tuple val(meta), val(pfsearch_out)
        val dirpath
        val blacklist_file

    output:
        tuple val(meta), path("prositeprofiles.json")

    exec:
    Map matches = [:]
    SignatureLibraryRelease library = new SignatureLibraryRelease("PROSITE profiles", null)
    def toSkip = new File("${dirpath.toString()}/${blacklist_file}").readLines()

    new File(pfsearch_out.toString()).eachLine { line ->
        if (line.trim()) {
            line = line.split()
            assert line.size() == 10
            String modelAccession = line[0].split("\\|")[0]
            if (modelAccession in toSkip) {
                return // skip flagged accessions
            }
            def matchData = getMatchData(line)
            seqId = matchData.sequenceId
            matches.computeIfAbsent(seqId) { [:] }
            Match matchObj = matches[seqId].computeIfAbsent(modelAccession) {
                new Match(matchData.profile, new Signature(modelAccession, library))
            }
            String alignment = matchData.alignment
            String cigarAlignment = Match.encodeCigarAlignment(alignment)
            Location location = new Location(matchData.start, matchData.end, matchData.normScore, alignment, cigarAlignment)
            matchObj.addLocation(location)
        }
    }
    def outputFilePath = task.workDir.resolve("prositeprofiles.json")
    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}

Map<String, Object> getMatchData(line) {
    def profile = line[0].split("\\|")
    return [
        profile     : profile[0],
        profileName : profile[1],
        sequenceId  : line[3],
        start       : line[4].toInteger(),
        end         : line[5].toInteger(),
        normScore   : line[7].toFloat(),
        alignment   : line[9]
    ]
}
