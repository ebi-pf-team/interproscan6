import groovy.json.JsonSlurper
import uk.ac.ebi.interpro.ProcessOutputGFF3
import uk.ac.ebi.interpro.ProcessOutputJSON
import uk.ac.ebi.interpro.ProcessOutputTSV
import uk.ac.ebi.interpro.ProcessOutputXML

if (args.size() != 7) {
    System.err.println "Usage: write-output.groovy <format> <input-directory> <metadata.json> <database.sqlite> <is-nucleic> <interproscan-version> <output.json>"
    System.exit(1)
}

def format = args[0]
def inputDir = new File(args[1])
def metadataPath = args[2]
def databasePath = args[3]
def isNucleic = args[4] == "true"
def iprscanVersion = args[5]
def outputPath = args[6]

def inputPaths = []
inputDir.eachFileRecurse { file ->
    if (file.name.endsWith('.json') && file.name != metadataPath) {
        inputPaths << file.absolutePath
    }
}

assert ["gff3", "json", "jsonl", "tsv", "xml"].contains(format)
assert inputPaths.size() > 0

def dbReleases = null
if (metadataPath != "-") {
    dbReleases = new JsonSlurper().parse(new File(metadataPath))

} else {
    dbReleases = null
}

if (format == "gff3") {
    ProcessOutputGFF3.run(inputPaths, databasePath, isNucleic, iprscanVersion, outputPath)
} else if (format == "json") {
    ProcessOutputJSON.run(inputPaths, databasePath, dbReleases, isNucleic, iprscanVersion, false, outputPath)
} else if (format == "jsonl") {
    ProcessOutputJSON.run(inputPaths, databasePath, dbReleases, isNucleic, iprscanVersion, true, outputPath)
} else if (format == "tsv") {
    ProcessOutputTSV.run(inputPaths, databasePath, isNucleic, outputPath)
} else if (format == "xml") {
    ProcessOutputXML.run(inputPaths, databasePath, dbReleases, isNucleic, iprscanVersion, outputPath)
}
