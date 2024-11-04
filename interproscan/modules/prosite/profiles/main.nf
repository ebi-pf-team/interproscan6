import groovy.io.FileType
import groovy.json.JsonOutput

process RUN_PFSEARCH {
    label 'prosite_pfsearch_runner'

    input:
        tuple val(meta), path(fasta)
        val models_dir

    output:
        tuple val(meta), path("prosite_profiles.out")

    script:
    """
    touch prosite_profiles.out
    find ${models_dir} -type f | while read profile; do
        output=\$(/opt/pftools/pfsearchV3 "\${profile}" "${fasta}" -f -o 7 -t 4)
        if [[ -n "\${output}" ]]; then
            echo "\${output}" >> prosite_profiles.out
        fi
    done
    """
}

process PARSE_PFSEARCH {
    label 'analysis_parser'

    input:
        tuple val(meta), val(pfsearch_out)
        val skip_flagged_profiles

    output:
        tuple val(meta), path("pfsearch_parsed.json")

    exec:
    Map matches = [:]
    def toSkip = new File(skip_flagged_profiles).readLines()

    new File(pfsearch_out.toString()).eachLine { line ->
        if (line.trim()) {
            line = line.split()
            assert line.size() == 10
            String modelAccession = line[0].split("\\|")[0]
            if (modelAccession in toSkip) {
                // skip flagged accessions
            }
            def match = createPrositeMatch(line)
            matchObj = matches.computeIfAbsent(match.sequenceId) { new Match(match.profile) }
            Location location = new Location(match.start, match.end, match.normScore, match.alignment)
            matchObj.addLocation(location)
        }
    }
    def outputFilePath = task.workDir.resolve("pfsearch_parsed.json")
    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}

def createPrositeMatch(line) {
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
