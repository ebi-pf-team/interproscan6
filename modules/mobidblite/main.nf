import groovy.json.JsonOutput

process RUN_MOBIDBLITE {
    label 'mem_min', 'time_short', 'dynamic', 'mobidblite_container'

    input:
    tuple val(meta), val(meta2), path(fasta)

    output:
    tuple val(meta), val(meta2), path("output.tsv")

    script:
    """
    idrpred --tempdir . --threads ${task.cpus} ${fasta} output.tsv
    """
}


process PARSE_MOBIDBLITE {
    label    'mem_low', 'time_short'
    executor 'local'

    input:
    tuple val(meta), val(meta2), val(mobidblite_output)

    output:
    tuple val(meta), val(meta2), path("mobidblite.json")

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
            SignatureLibraryRelease library = new SignatureLibraryRelease("MobiDB-lite", "4.0")
            Signature signature = new Signature("mobidb-lite", "disorder_prediction", "consensus disorder prediction", library, null)
            match = new Match("mobidb-lite", signature)
            matches[sequenceId] = [:]
            matches[sequenceId]["mobidb-lite"] = match
        }

        match.addLocation(new Location(start, end, feature))
    }

    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}
