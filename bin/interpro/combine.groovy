import uk.ac.ebi.interpro.Process

if (args.size() != 2) {
    System.err.println "Usage: combine.groovy <input_dir> <output.json>"
    System.exit(1)
}

def inputDir = new File(args[0])
def outputPath = args[1]
def inputPaths = []
inputDir.eachFileRecurse { file ->
    if (file.name.endsWith('.json')) {
        inputPaths << file.absolutePath
    }
}

assert inputPaths.size() > 0
Process.combineMatches(inputPaths, outputPath)
