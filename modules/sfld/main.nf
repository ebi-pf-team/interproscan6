import groovy.json.JsonOutput

process SEARCH_SFLD {
    label 'small', 'ips6_container'

    input:
    tuple val(meta), path(fasta)
    path dirpath
    val hmmfile

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
        ${dirpath}/${hmmfile} ${fasta}
    """
}

process POST_PROCESS_SFLD {
    label 'small'

    input:
    tuple val(meta), path(hmmsearch_out), val(hmmsearch_dtbl), val(hmmsearch_alignment)
    path dirpath
    val annofile

    output:
    tuple val(meta), path("sfld.tsv"), path(hmmsearch_out)

    script:
    """
    ${projectDir}/bin/sfld/sfld_postprocess \
        --alignment "${hmmsearch_alignment}" \
        --dom "${hmmsearch_dtbl}" \
        --hmmer-out "${hmmsearch_out}" \
        --site-info "${dirpath}/${annofile}" \
        --output sfld.tsv
    # ${hmmsearch_out} "hmmsearch.out"
    """
}

process PARSE_SFLD {
    label 'native'

    input:
    tuple val(meta), val(postprocess_out), val(hmmsearch_out)
    val dirpath
    val hierarchydb

    output:
    tuple val(meta), path("sfld.json")

    exec:
    def outputFilePath = task.workDir.resolve("sfld.json")
    def (hmmLengths, hmmBounds) = getHmmData(hmmsearch_out.toString())
    def sequences = parseOutput(postprocess_out.toString(), hmmLengths, hmmBounds)
    def hierarchy = parseHierarchy("${dirpath.toString()}/${hierarchydb}")
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
                selectedMatches.addAll(promotedMatches)
            }
        }
        
        if (selectedMatches.size() > 0) {
            /*
                Cases where matches have the same locations and the same ancestor (e.g., SFLDF00273
                and SFLDF00413 have the same ali location and the same ancestor — SFLDS00029 —
                so SFLDS00029 is duplicated when promoted).
            */
            List<Match> uniqueMatches = []
            Set<String> seenKeys = [] as Set
            // sorting matches by location evalue ASC, location score DESC to keep the best matches
            List<Match> sortedMatches = selectedMatches.sort { a, b ->
                (a.locations[0].evalue <=> b.locations[0].evalue) ?: -(a.locations[0].score <=> b.locations[0].score)
            }
            sortedMatches.each { match ->
                String key = "${match.modelAccession}:${match.locations[0].start}:${match.locations[0].end}"
                if (seenKeys.contains(key)) {
                    return
                }
                uniqueMatches.add(match)
                seenKeys.add(key)
            }
            selectedMatches = uniqueMatches
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
            def childAccession = accessions[-2]  // the final component is the description, only take the accs
            def parents = accessions.subList(0, accessions.size() - 2)
            hierarchy[childAccession] = parents as Set
            for (parent in parents) {  // ensure all ancestors are in the hierarchy
                ancestors = hierarchy.get(parent)
                if (ancestors.size() > 0) {
                    hierarchy[childAccession].addAll(ancestors)
                }
            }
	    }
    }
    return hierarchy
}

List<Map> getHmmData(String outputFilePath) {
    def matchesMap = HMMER3.parseOutput(outputFilePath.toString(), "AntiFam")  // proteinMd5: modelAccession: Match
    def hmmLengths = [:]  // modelAccession: hmmLength
    def hmmBounds  = [:]  // proteinMd5: modelAccession: (tuple representing loc): hmmbound
    matchesMap.each { String sequenceId, matches ->
        matches.each { String modelAccession, Match match ->
            match.locations.each { Location loc ->
                def locKey = [loc.start, loc.end, loc.hmmStart, loc.hmmEnd, loc.envelopeStart, loc.envelopeEnd]
                hmmLengths[modelAccession] = loc.hmmLength
                hmmBounds.computeIfAbsent(sequenceId) {     [:] }
                         .computeIfAbsent(modelAccession) { [:] }
                         .computeIfAbsent(locKey) {         loc.hmmBounds }
            }
        }
    }
    return [hmmLengths, hmmBounds]
}

Map<String, Map<String, Match>> parseOutput(
    String outputFilePath,
    Map<String, Integer> hmmLengthsMap,
    Map<String, Map> hmmBoundsMap
) {
    // Parse the output TSV file from the SFLD postprocess bin
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

            matches[sequenceId] = parseBlock(reader, sequenceId, hmmLengthsMap, hmmBoundsMap)
        }
    }

    return matches
}

Map<String, Match> parseBlock(
    Reader reader,
    String sequenceId,
    Map<String, Integer> hmmLengthsMap,
    Map<String, Map> hmmBoundsMap     // proteinMd5: modelAccession: (tuple representing loc): hmmbound
) {
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
            Integer hmmLength = hmmLengthsMap[modelAccession]
            def locKey = [aliStart, aliEnd, hmmStart, hmmEnd, envStart, envEnd]
            String hmmBounds = (hmmBoundsMap.containsKey(sequenceId) &&
                    hmmBoundsMap[sequenceId].containsKey(modelAccession) &&
                    hmmBoundsMap[sequenceId][modelAccession] != null) ?
                    hmmBoundsMap[sequenceId][modelAccession][locKey]?.toString() :
                    null

            if (aliStart <= aliEnd && hmmStart <= hmmEnd && envStart <= envEnd) {
                Match match = domains.computeIfAbsent(modelAccession, k -> {
                    Signature signature = new Signature(modelAccession, library)
                    return new Match(modelAccession, seqEvalue, seqScore, seqBias, signature)
                });
                Location location = new Location(aliStart, aliEnd, hmmStart, hmmEnd, hmmLength, hmmBounds,
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
