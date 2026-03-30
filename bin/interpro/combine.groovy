import java.nio.file.Path
import java.nio.file.Paths
import uk.ac.ebi.interpro.ProcessCombine

if (args.size() != 2) {
    System.err.println "Usage: combine.groovy <input_dir> <output.json>"
    System.exit(1)
}

def inputDir = Paths.get(args[0])
def outputPath = Paths.get(args[1])
def inputPaths = []
inputDir.eachFileRecurse { file ->
    if (file.fileName.toString().endsWith('.json')) {
        inputPaths << file
    }
}

assert inputPaths.size() > 0
ProcessCombine.run(inputPaths, outputPath)
