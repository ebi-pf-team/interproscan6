#!/bin/bash
//bin/true; CLASSPATH=$1; shift; exec groovy -cp "$CLASSPATH:$CLASSPATH/*:." "$0" "$@"

import java.nio.file.Path
import java.nio.file.Paths
import uk.ac.ebi.interpro.ProcessReprLocations

if (args.size() != 2) {
    System.err.println "Usage: select-repr-locations.groovy <input-dir> <output-dir>"
    System.exit(1)
}

Path inputDir = Paths.get(args[0])
Path outputDir = Paths.get(args[1])
def inputPaths = []
inputDir.eachFileRecurse { file ->
    if (file.fileName.toString().endsWith('.json')) {
        inputPaths << file
    }
}

ProcessReprLocations.run(inputPaths, outputDir)
