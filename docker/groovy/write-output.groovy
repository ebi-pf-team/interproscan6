import groovy.json.JsonSlurper
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import uk.ac.ebi.interpro.ProcessOutputGFF3
import uk.ac.ebi.interpro.ProcessOutputJSON
import uk.ac.ebi.interpro.ProcessOutputTSV
import uk.ac.ebi.interpro.ProcessOutputXML

if (args.size() != 7) {
    System.err.println "Usage: write-output.groovy <format> <input-directory> <database.sqlite> <is-nucleic> <interproscan-version> <interpro-version> <output-file>"
    System.exit(1)
}

def format = args[0]
def inputDir = Paths.get(args[1])
def databasePath = Paths.get(args[2])
def isNucleic = args[3] == "true"
def iprscanVersion = args[4]
def interproVersion = args[5]
def outputPath = Paths.get(args[6])
def outputPathNormalized = outputPath.toAbsolutePath().normalize()

def inputPaths = [] 
inputDir.eachFileRecurse { file ->
    if (file.fileName.toString().endsWith('.json') &&
        !file.toAbsolutePath().normalize().equals(outputPathNormalized)) {
        inputPaths << file
    }
}

assert ["gff3", "json", "jsonl", "tsv", "xml"].contains(format)
assert inputPaths.size() > 0

if (format == "gff3") {
    ProcessOutputGFF3.run(inputPaths, databasePath, isNucleic, iprscanVersion, outputPath)
} else if (format == "json") {
    ProcessOutputJSON.run(inputPaths, databasePath, isNucleic, iprscanVersion, interproVersion, false, outputPath)
} else if (format == "jsonl") {
    ProcessOutputJSON.run(inputPaths, databasePath, isNucleic, iprscanVersion, interproVersion, true, outputPath)
} else if (format == "tsv") {
    ProcessOutputTSV.run(inputPaths, databasePath, isNucleic, outputPath)
} else if (format == "xml") {
    ProcessOutputXML.run(inputPaths, databasePath, isNucleic, iprscanVersion, interproVersion, outputPath)
}
