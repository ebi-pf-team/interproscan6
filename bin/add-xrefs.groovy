#!/bin/bash
//bin/true; CLASSPATH=$1; shift; exec groovy -cp "$CLASSPATH:$CLASSPATH/*:." "$0" "$@"

import java.nio.file.Path
import java.nio.file.Paths
import uk.ac.ebi.interpro.ProcessXrefs

if (args.size() != 6) {
    System.err.println "Usage: add-xrefs.groovy <input-dir> <interpro-dir> <panther-paint-dir> <add-go-terms> <add-pathways> <output-dir>"
    System.exit(1)
}

Path inputDir = Paths.get(args[0])
Path interproDir = args[1] == "-" ? null : Paths.get(args[1])
Path pantherPaintDir = args[2] == "-" ? null : Paths.get(args[2])
def addGoTerms = args[3] == "true"
def addPathways = args[4] == "true"
Path outputDir = Paths.get(args[5])

def inputPaths = []
inputDir.eachFileRecurse { file ->
    if (file.fileName.toString().endsWith('.json')) {
        inputPaths << file
    }
}

ProcessXrefs.run(inputPaths, interproDir, pantherPaintDir, addGoTerms, addPathways, outputDir)
