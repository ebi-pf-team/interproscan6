import groovy.json.JsonSlurper
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import uk.ac.ebi.interpro.ProcessOutputGFF3
import uk.ac.ebi.interpro.ProcessOutputJSON
import uk.ac.ebi.interpro.ProcessOutputTSV
import uk.ac.ebi.interpro.ProcessOutputXML

if (args.size() != 7) {
    System.err.println "Usage: write-output.groovy <format> <input-directory> <metadata.json> <database.sqlite> <is-nucleic> <interproscan-version> <output.json>"
    System.exit(1)
}

def format = args[0]
Path inputDir = Paths.get(args[1])
Path metadataPath = args[2] == "-" ? null : Paths.get(args[2])
Path databasePath = Paths.get(args[3])
def isNucleic = args[4] == "true"
def iprscanVersion = args[5]
Path outputPath = Paths.get(args[6])

List<Path> inputPaths = []
Path normalizedMetadataPath = metadataPath?.toAbsolutePath()?.normalize()
Files.walk(inputDir).withCloseable { inputDirPaths -> 
    inputDirPaths.forEach { Path file ->
        if (Files.isRegularFile(file) &&
            file.fileName.toString().endsWith('.json') &&
            file.toAbsolutePath().normalize() != normalizedMetadataPath) {
            inputPaths << file
        }
    }
}

assert ["gff3", "json", "jsonl", "tsv", "xml"].contains(format)
assert inputPaths.size() > 0

def dbReleases = null
if (metadataPath != null) {
    dbReleases = metadataPath.newReader().withCloseable { reader -> 
        new JsonSlurper().parse(reader)
    }
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
