import java.nio.file.Path
import java.nio.file.Paths
import groovy.json.JsonSlurper
import uk.ac.ebi.interpro.ProcessXrefs

if (args.size() != 6) {
    System.err.println "Usage: add-xrefs.groovy <input.json> <interpro-dir> <panther-paint-dir> <add-go-terms> <add-pathways> <output.json>"
    System.exit(1)
}

Path inputPath = Paths.get(args[0])
Path interproDir = args[1] == "-" ? null : Paths.get(args[1])
Path pantherPaintDir = args[2] == "-" ? null : Paths.get(args[2])
def addGoTerms = args[3] == "true"
def addPathways = args[4] == "true"
Path outputPath = Paths.get(args[5])

ProcessXrefs.run(inputPath, interproDir, pantherPaintDir, addGoTerms, addPathways, outputPath)
