import groovy.json.JsonOutput

process RUN_SFLD {
    label 'small', 'ips6_container'

    input:
    tuple val(meta), path(fasta)
    path hmmdb

    output:
    tuple val(meta), path("hmmsearch.out"), path("hmmsearch.tab"), path("hmmsearch.sto")

    script:
    """
    /opt/hmmer3/bin/hmmsearch \
        -Z 378 --acc \
        --cut_ga \
        --cpu ${task.cpus} \
        -o hmmsearch.out \
        --domtblout hmmsearch.tab \
        -A hmmsearch.sto \
        ${hmmdb} ${fasta}
    """
}

process POST_PROCESS_SFLD {
    label 'small'

    input:
    tuple val(meta), val(hmmsearch_out), val(hmmsearch_dtbl), val(hmmsearch_alignment)
    val site_info

    output:
    tuple val(meta), path("sfld.tsv")

    script:
    """
    $projectDir/interproscan/bin/sfld/sfld_postprocess \
        --alignment ${hmmsearch_alignment} \
        --dom ${hmmsearch_dtbl} \
        --hmmer-out ${hmmsearch_out} \
        --site-info ${site_info} \
        --output sfld.tsv
    """
}


process PARSE_SFLD {
    label 'local'

    input:
    tuple val(meta), val(postprocess_out)
    val hierarchy_db

    output:
    tuple val(meta), path("sfld.json")

    exec:
    def outputFilePath = task.workDir.resolve("sfld.json")
    def sequences = parseOutput(postprocess_out.toString())
    def hierarchy = parseHierarchy(hierarchy_db.toString())
    SignatureLibraryRelease library = new SignatureLibraryRelease("SFLD", null)

    sequences = sequences.collectEntries { seqId, matches -> 
        // Flatten matches (one location per match)
        matches = matches.collectMany { key, match ->
            return match.locations.collect { location ->
                Signature signature = new Signature(match.modelAccession, library)
                Match newMatch = new Match(match.modelAccession, match.evalue, match.score, match.bias, signature)
                newMatch.addLocation(location.clone())
                return newMatch
            }
        }

        if (matches.size() > 1) {
            // Now resolve matches to remove overlapping matches from the same hierarchy
            Set<Integer> ignored = [] as Set
            matches = matches
                .collect { match ->
                    if (match.modelAccession.startsWith("SFLDF")) {
                        // Hihgly specific model: always add
                        return match
                    }

                    boolean overlaps = false
                    for (Match otherMatch: matches) {
                        if (ignored.contains(otherMatch)) {
                            // Current match overlaps with another match: ignore
                            continue
                        } else if (match.modelAccession == otherMatch.modelAccession) {
                            // Two matches/locations from the same model: keep
                            continue
                        }

                        Location l1 = match.locations[0]
                        Location l2 = otherMatch.locations[0]

                        if (!(l1.start > l2.end || l2.start > l1.end)) {
                            // Matches overlap
                            
                            def otherParents = hierarchy.get(otherMatch.modelAccession)
                            if (otherParents != null && otherParents.contains(match.modelAccession)) {
                                // Current match overlaps a more specific match: we want to keep the most specific
                                overlaps = true
                                break
                            }    
                        }
                    }

                    if (overlaps) {
                        // Mark this match as to be ignored
                        ignored.add(match)
                        return null
                    }

                    return match
                }
                .findAll { it != null }
        }

        def selectedMatches = [] as Set
        matches.each { match ->
            def parents = hierarchy.get(match.modelAccession)
            if (parents) {
                /*
                    Propagate match up through it hierarchy.
                    If there are Superfamily, Group, and Family models in a tree,
                    and a sequence matches F, it should inherit the S and G annotations
                */
                def promotedMatches = parents
                    .findAll { it != match.modelAccession }
                    .collect {
                        Signature signature = new Signature(it, library)
                        Match promotedMatch = new Match(it, match.evalue, match.score, match.bias, signature)
                        promotedMatch.addLocation(match.locations[0].clone())
                        return promotedMatch
                    }

                /*
                    Check newly promoted matches against previously selected matches:
                        - promoted fully contains other match: keep promoted and remove other
                        - promoted fully within other match: skip promoted
                        - partial or no overlap: keep promoted
                */
                boolean keepPromoted = true
                Set<Match> toRemove = [] as Set
                promotedMatches.each { promotedMatch ->
                    for (Match selectedMatch: selectedMatches) {
                        if (promotedMatch.locations[0].start <= selectedMatch.locations[0].start &&
                            promotedMatch.locations[0].end >= selectedMatch.locations[0].end) {
                            toRemove.add(selectedMatch)
                        } else if (promotedMatch.locations[0].start >= selectedMatch.locations[0].start &&
                                  promotedMatch.locations[0].end <= selectedMatch.locations[0].end) {
                            keepPromoted = false
                            break
                        }
                    }

                    toRemove.each { selectedMatches.remove(it) }

                    if (keepPromoted) {
                        selectedMatches.add(promotedMatch)
                    }
                }
            }
        }
        
        if (selectedMatches.size() > 0) {
            /*
            KEEP THIS BLOCK
            I5 has a "remove duplicates" step, but it seems bugged,
            and I couldn't even find one case where it had consequence.
            We'll implement the same thing here if we want such a case.
            */
        }

        // Add initial matches (the ones used for promotion)
        selectedMatches.addAll(matches)

        // List to Map
        def finalMatches = [:]
        selectedMatches.each { match ->
            Match finalMatch = finalMatches.get(match.modelAccession)
            if (finalMatch) {
                finalMatch.addLocation(match.locations[0])
            } else {
                finalMatches[match.modelAccession] = match
            }
        }

        return [seqId, finalMatches]
    }

    def json = JsonOutput.toJson(sequences)
    new File(outputFilePath.toString()).write(json)
}

Map<String, Set<String>> parseHierarchy(String filePath) {
    def hierarchy = [:]
    File file = new File(filePath)
    file.withReader{ reader ->
        for (String line = reader.readLine(); line != null; line = reader.readLine()) {
            def accessions = line.split(/\t/).toList()
            def childAccession = accessions[-2]  // the last one is the description
            def parents = accessions.subList(0, accessions.size() - 2)
            hierarchy[childAccession] = parents as Set
            for (parent in parents) {  // making sure all the ancestors are in the hierarchy
                ancestors = hierarchy.get(parent)
                if (ancestors.size() > 0) {
                    hierarchy[childAccession].addAll(ancestors)
                }
            }
	}
    }
    return hierarchy
}

Map<String, Map<String, Match>> parseOutput(String outputFilePath) {
    def matches = [:]
    File file = new File(outputFilePath)
    file.withReader{ reader ->
        while (true) {
            String sequenceId = null

            for (String line = reader.readLine(); line != null; line = reader.readLine()) {
                def matcher = line.trim() =~ ~/^Sequence:\s+(\S+)$/
                if (matcher.find()) {
                    sequenceId = matcher.group(1).trim()
                    break
                }
            }

            if (sequenceId == null) {
                break
            }

            matches[sequenceId] = parseBlock(reader)
        }
    }

    return matches
}

Map<String, Match> parseBlock(Reader reader) {
    SignatureLibraryRelease library = new SignatureLibraryRelease("SFLD", null)
    boolean inDomains = false
    def domains = [:]
    while (true) {
        String line = reader.readLine()?.trim()
        if (!line) break
        if (line == "Domains:") {
            inDomains = true
        } else if (line == "Sites:") {
            inDomains = false
        } else if (line == "//") {
            break
        } else if (inDomains) {
            def fields = line.split(/\t/)
            assert fields.length == 15
            String modelAccession = fields[0]
            Double seqEvalue = fields[1] as Double
            Double seqScore = fields[2] as Double
            Double seqBias = fields[3] as Double
            Integer hmmStart = fields[4] as Integer
            Integer hmmEnd = fields[5] as Integer
            Double domScore = fields[6] as Double
            int aliStart = fields[7] as int
            int aliEnd = fields[8] as int
            Integer envStart = fields[9] as Integer
            Integer envEnd = fields[10] as Integer
            Double domCEvalue = fields[11] as Double
            Double domIEvalue = fields[12] as Double
            Double accuracy = fields[13] as Double
            Double domBias = fields[14] as Double

            if (aliStart <= aliEnd && hmmStart <= hmmEnd && envStart <= envEnd) {
                Match match = domains.computeIfAbsent(modelAccession, k -> {
                    Signature signature = new Signature(modelAccession, library)
                    return new Match(modelAccession, seqEvalue, seqScore, seqBias, signature)
                });
                Location location = new Location(aliStart, aliEnd, hmmStart, hmmEnd, null, null,
                        envStart, envEnd, domIEvalue, domScore, domBias)
                match.addLocation(location)
            }
        } else {
            def fields = line.split(/\s/, 3)
            assert fields.length == 2 || fields.length == 3
            String modelAccession = fields[0]
            String residues = fields[1]
            String description = fields.length == 3 ? fields[2] : null

            Match match = domains.get(modelAccession)
            if (match != null) {
                Site site = new Site(description, residues)
                match.addSite(site)                        
            }
        }
    }

    return domains
}
