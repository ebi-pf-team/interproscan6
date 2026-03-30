import java.nio.file.Paths
import uk.ac.ebi.interpro.ProcessReprLocations

if (args.size() != 2) {
    System.err.println "Usage: select-repr-locations.groovy <input.json> <output.json>"
    System.exit(1)
}

ProcessReprLocations.run(Paths.get(args[0]), Paths.get(args[1]))
