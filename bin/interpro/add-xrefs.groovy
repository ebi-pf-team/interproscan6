import groovy.json.JsonSlurper
import uk.ac.ebi.interpro.Process

if (args.size() != 6) {
    System.err.println "Usage: add-xrefs.groovy <input.json> <metadata.json> <add-go-terms> <add-pathways> <panther-paint-dir> <output.json>"
    System.exit(1)
}

def inputPath = args[0]
def metadataPath = args[1]
def addGoTerms = args[2] == "true"
def addPathways = args[3] == "true"
def pantherPaintDirectory = args[4] 
def outputPath = args[5]

def dbReleases = new JsonSlurper().parse(new File(metadataPath))

Process.addXrefs(inputPath, dbReleases, addGoTerms, addPathways, pantherPaintDirectory, outputPath)
