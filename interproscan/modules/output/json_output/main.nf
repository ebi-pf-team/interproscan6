import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import java.util.regex.Pattern

// NT_PATTERN = Pattern.compile(/^orf\d+\s+source=(.*)\s+coords=.*$/)

process JSON_OUTPUT {
    label 'write_output'

    input:
    val seqMatches
    val outputPath
    val ips6Version

    exec:
    def jsonSlurper = new JsonSlurper()
    def jsonOutput = [:]
    jsonOutput["interproscan-version"] = ips6Version
    jsonOutput["results"] = []

    def matches = jsonSlurper.parse(seqMatches).each { sequence ->
        jsonOutput["results"].add(sequence)
    }

    def outputFilePath = "${outputPath}.ips6.json"
    def json = JsonOutput.toJson(jsonOutput)
    new File(outputFilePath.toString()).write(json)
}
