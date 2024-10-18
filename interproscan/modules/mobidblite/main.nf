import groovy.json.JsonOutput

process RUN_MOBIDBLITE {
    label 'mobidblite_runner'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("output.tsv")

    script:
    """
    idrpred ${fasta} output.tsv
    """
}


process PARSE_MOBIDBLITE {
    label 'analysis_parser'

    input:
    tuple val(meta), val(mobidblite_output)

    output:
    tuple val(meta), path("mobidblite.json")

    exec:
    def outputFilePath = task.workDir.resolve("mobidblite.json")
    Match match = null
    def matches = [:]
    file(mobidblite_output.toString()).eachLine { line ->
        def lineData = line.split(/\s+/)
        def sequenceId = lineData[0]
        def start = lineData[1].toInteger()
        def end = lineData[2].toInteger()
        def feature = lineData[3] != "-" ? lineData[3] : null

        if (matches.containsKey(sequenceId)) {
            match = matches[sequenceId]["mobidb-lite"]
        } else {
            match = new Match("mobidb_lite")
            matches[sequenceId] = [:]
            matches[sequenceId]["mobidb-lite"] = match
        }

        match.addLocation(new Location(start, end, feature))
    }

    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}
