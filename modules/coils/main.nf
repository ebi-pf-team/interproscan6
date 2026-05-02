process RUN_COILS {
    label     'mem_min', 'time_veryshort'
    container 'interpro/ncoils:2.2.1'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("coils.out")

    script:
    """
    ncoils -c < ${fasta} > coils.out
    """
}


process PARSE_COILS {
    label    'mem_low', 'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(coils_out)

    output:
    tuple val(meta), path("coils.json")

    exec:
    def matches = [:]
    def sequenceId = null
    def library = new uk.ac.ebi.interpro.SignatureLibraryRelease("COILS", "2.2.1")

    coils_out.eachLine { line ->
        line = line.trim()
        if (line.startsWith(">")) {
            // Coils report the full sequence header (ID + description)
            sequenceId = line.substring(1).split()[0]
            matches[sequenceId] = [:]

            def match = new uk.ac.ebi.interpro.Match("Coil", new uk.ac.ebi.interpro.Signature("Coil", "Coil", null, library, null))
            matches[sequenceId]["Coil"] = match
        } else if (line != "//" && sequenceId) {
            def fields = line.split(/\s+/)
            def start = fields[0].toInteger()
            def end = fields[1].toInteger()
            matches[sequenceId]["Coil"].addLocation(new uk.ac.ebi.interpro.Location  (start, end))
        }
    }

    def filepath = task.workDir.resolve("coils.json")
    filepath.text = groovy.json.JsonOutput.toJson(
        matches.findAll { m -> m.value["Coil"].locations.size() > 0 }
    )
}
