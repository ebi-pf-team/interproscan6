import uk.ac.ebi.interpro.ProcessReprLocations

if (args.size() != 2) {
    System.err.println "Usage: select-repr-locations.groovy <input.json> <output.json>"
    System.exit(1)
}

ProcessReprLocations.run(args[0], args[1])
