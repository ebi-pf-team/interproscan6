import groovy.io.FileType
import groovy.json.JsonOutput

process PFSEARCH_RUNNER {
    label 'prosite_pfsearch_runner'

    input:
        tuple val(meta), path(fasta)
        val models_dir
        val skip_flagged_profiles

    output:
        tuple val(meta), path("prosite_profiles.out")
        val skip_flagged_profiles

    exec:
    // Run pfsearchV3 for all PROSITE profiles
    def outputFilePath = task.workDir.resolve("prosite_profiles.out")
    def outputFile = new File(outputFilePath.toString())
    outputFile.createNewFile() // create the file to ensure an output file even if empty

    List<String> profilePaths = []
    new File(models_dir).eachFileRecurse(FileType.FILES) { file ->
        profilePaths << file.absolutePath
    }
    profilePaths.each { profile ->
        def runCmd = "/bin/bash -c '/opt/pftools/pfsearchV3 ${profile} ${fasta} -f -o 7 -t 4'"
        try {
            def process = runCmd.execute()
            def output = process.text
            if (output.trim()) {
               outputFile.append(output + '\n')
            }
        } catch (Exception e) {
            println "Error running pfsearchV3 using cmd: ${runCmd}\nError: ${e.message}"
            System.exit(1)
        }
    }
}

process PFSEARCH_PARSER {
    label 'analysis_parser'

    input:
        tuple val(meta), val(pfsearch_out)
        val skip_flagged_profiles

    output:
        tuple val(meta), path("pfsearch_parsed.json")

    exec:
    Map profilesMatches = [:]
    def toSkip = new File(skip_flagged_profiles).readLines()
    def DATA_LINE = ~/^([\w\|]+)\s+(\d+)\s+(\-?\d+)\s+([\d\)\(\w+\|\.\-\_]+|\s+)\s+(\d+)\s+(\d+)\s+(\-?\d+)\s+(\d+\.?\d+)\s+([\+\-])\s+([\w\-]+)$/

    new File(pfsearch_out.toString()).eachLine { line ->
        if (line.trim()) {
            def lineData = DATA_LINE.matcher(line)

            if (!lineData.matches()) {
                println "ERROR: the following line did not match the expected architecture"
                println line
            } else {
                def profileId = lineData.group(1).split("\\|")[0]
                if (profileId in toSkip) {
                    return
                }
                def hit = createPrositeHit(lineData)

                if (patternsMatches.containsKey(hit.sequenceId)) {
                    match = profilesMatches[hit.sequenceId]
                } else {
                    match = new Match(hit.profile)
                    match.name = hit.profileName
                    profilesMatches[hit.sequenceId] = match
                }

                Location location = new Location(hit.start, hit.end, hit.normScore, hit.alignment)
                match.addLocation(location)
            }
        }
    }

    def outputFilePath = task.workDir.resolve("pfsearch_parsed.json")
    def json = JsonOutput.toJson(profilesMatches)
    new File(outputFilePath.toString()).write(json)
}

def createPrositeHit(value) {
    def profile, profileName = value[0].split("\\|")
    return [
        profile     : profile,
        profileName : profileName,
        motifStart  : value[1].toInteger(),
        motifEnd    : value[2].toInteger(),
        sequenceId  : value[3],
        start       : value[4].toInteger(),
        end         : value[5].toInteger(),
        rawScore    : value[6].toFloat(),
        normScore   : value[7].toFloat(),
        symbol      : value[8],
        alignment   : value[9]
    ]
}
