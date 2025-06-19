import groovy.json.JsonOutput

import Match

process RUN_COILS {
    label 'mini', 'ips6_container'

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
    label    'tiny'
    executor 'local'

    input:
    tuple val(meta), val(coils_out)

    output:
    tuple val(meta), path("coils.json")

    exec:
    def outputFilePath = task.workDir.resolve("coils.json")
    def matches = [:]
    def sequenceId = null
    SignatureLibraryRelease library = new SignatureLibraryRelease("COILS", "2.2.1")

    file(coils_out.toString()).eachLine { line ->
        line = line.trim()
        if (line.startsWith(">")) {
            // Coils report the full sequence header (ID + description)
            sequenceId = line.substring(1).split()[0]
            matches[sequenceId] = [:]

            Match match = new Match("Coil", new Signature("Coil", "Coil", null, library, null))
            matches[sequenceId]["Coil"] = match
        } else if (line != "//" && sequenceId) {
            def fields = line.split(/\s+/)
            def start = fields[0].toInteger()
            def end = fields[1].toInteger()
            matches[sequenceId]["Coil"].addLocation(new Location(start, end))
        }
    }

    def json = JsonOutput.toJson(matches.findAll { it.value["Coil"].locations.size() > 0 })
    new File(outputFilePath.toString()).write(json)
}
