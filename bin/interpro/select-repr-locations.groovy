import uk.ac.ebi.interpro.Process

if (args.size() != 2) {
    System.err.println "Usage: select-repr-locations.groovy <input.json> <output.json>"
    System.exit(1)
}

Process.selectRepresentativeLocations(args[0], args[1])
