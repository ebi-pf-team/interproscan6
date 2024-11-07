import groovy.json.JsonOutput

process RUN_COILS {
    label 'coils_runner'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("coils.out")

    script:
    """
    /opt/coils/ncoils -c < ${fasta} > coils.out
    """
}


process PARSE_COILS {
    label 'analysis_parser'

    input:
    tuple val(meta), val(coils_out)

    output:
    tuple val(meta), path("coils.json")

    exec:
    def outputFilePath = task.workDir.resolve("coils.json")
    def matches = [:]
    def sequenceId = null

    file(coils_out.toString()).eachLine { line ->
        line = line.trim()
        if (line.startsWith(">")) {
            // Coils report the full sequence header (ID + description)
            sequenceId = line.substring(1).split()[0]
            matches[sequenceId] = [:]
            matches[sequenceId]["Coil"] = new Match("coils")
        } else if (line != "//" && sequenceId) {
            def fields = line.split(/\s+/)
            def start = fields[0].toInteger()
            def end = fields[1].toInteger()
            matches[sequenceId]["Coil"].addLocation(new Location(start, end))
        }
    }

    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}
