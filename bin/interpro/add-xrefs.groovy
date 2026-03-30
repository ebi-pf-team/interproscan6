import java.nio.file.Path
import java.nio.file.Paths
import groovy.json.JsonSlurper
import uk.ac.ebi.interpro.ProcessXrefs

if (args.size() != 6) {
    System.err.println "Usage: add-xrefs.groovy <input.json> <metadata.json> <add-go-terms> <add-pathways> <panther-paint-dir> <output.json>"
    System.exit(1)
}

Path inputPath = Paths.get(args[0])
Path metadataPath = Paths.get(args[1])
def addGoTerms = args[2] == "true"
def addPathways = args[3] == "true"
Path pantherPaintDirectory = Paths.get(args[4] )
Path outputPath = Paths.get(args[5])

def dbReleases = new JsonSlurper().parse(metadataPath)

ProcessXrefs.run(inputPath, dbReleases, addGoTerms, addPathways, pantherPaintDirectory, outputPath)
