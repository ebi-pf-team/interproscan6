import groovy.json.JsonOutput

process RUN_MOBIDBLITE {
    label 'mini', 'dynamic', 'mobidblite_container'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("output.tsv")

    script:
    """
    idrpred --tempdir . --threads ${task.cpus} ${fasta} output.tsv
    """
}


process PARSE_MOBIDBLITE {
    label    'tiny'
    executor 'local'

    input:
    tuple val(meta), val(mobidblite_output)

    output:
    tuple val(meta), path("mobidblite.json")

    exec:
    def outputFilePath = task.workDir.resolve("mobidblite.json")
    Match match = null
    def matches = [:]
    file(mobidblite_output.toString()).eachLine { line ->
        def fields = line.split(/\t/)
        assert fields.size() == 4
        def sequenceId = fields[0]
        def start = fields[1].toInteger()
        def end = fields[2].toInteger()
        def feature = fields[3] != "-" ? fields[3] : null

        if (matches.containsKey(sequenceId)) {
            match = matches[sequenceId]["mobidb-lite"]
        } else {
            match = new Match("mobidb-lite")
            SignatureLibraryRelease library = new SignatureLibraryRelease("MobiDB-lite", "4.0")
            match.signature = new Signature("mobidb-lite", "disorder_prediction", "consensus disorder prediction", library, null)
            matches[sequenceId] = [:]
            matches[sequenceId]["mobidb-lite"] = match
        }

        match.addLocation(new Location(start, end, feature))
    }

    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}
