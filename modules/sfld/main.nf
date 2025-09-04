import groovy.json.JsonOutput

import Match

process SEARCH_SFLD {
    label 'mini', 'dynamic', 'ips6_container'

    input:
    tuple val(meta), path(fasta)
    path dirpath
    val hmmfile

    output:
    tuple val(meta), path("hmmsearch.out"), path("hmmsearch.tab"), path("hmmsearch.sto")

    script:
    """
    hmmsearch \
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
    label 'mini', 'ips6_container'

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
    label    'tiny'
    executor 'local'

    input:
    tuple val(meta), val(postprocess_out), val(hmmsearch_out)
    val dirpath
    val hierarchydb

    output:
    tuple val(meta), path("sfld.json")

    exec:
    def outputFilePath = task.workDir.resolve("sfld.json")
    def hmmerMatches = HMMER3.parseOutput(hmmsearch_out.toString(), "SFLD")
    def sequences = parseOutput(postprocess_out.toString(), hmmerMatches)
    def hierarchies = getHierarchies("${dirpath.toString()}/${hierarchydb}")
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
            // Remove overlapping matches from models in the same hierarchy (keep the most specific)
            Set<Integer> ignored = [] as Set
            matches = matches
                .collect { match ->
                    if (match.modelAccession.startsWith("SFLDF")) {
                        // Highly specific model (family): always keep
                        return match
                    }

                    boolean ignore = false
                    for (Match otherMatch: matches) {
                        if (match.modelAccession == otherMatch.modelAccession) {
                            // Two matches/locations from the same mode: keep
                            continue
                        } else if (ignored.contains(otherMatch)) {
                            // otherMatch has previously been flagged as to be ignored
                            continue
                        }

                        Location l1 = match.locations[0]
                        Location l2 = otherMatch.locations[0]

                        if (!(l1.start > l2.end || l2.start > l1.end)) {
                            // Matches overlap
                            
                            def otherParents = hierarchies.get(otherMatch.modelAccession)
                            if (otherParents != null && otherParents.contains(match.modelAccession)) {
                                /*
                                    The current match overlaps a match from a more specific model:
                                    we want to keep the most specific match, so we ignore the current match
                                */
                                ignore = true
                                break
                            }    
                        }
                    }

                    if (ignore) {
                        // Mark this match as to be ignored
                        ignored.add(match)
                        return null
                    }

                    return match
                }
                .findAll { it != null }
        }

        def matchesWithPromoted = matches.clone()
        matches.each { match ->
            def parents = hierarchies.get(match.modelAccession)
            if (parents) {
                /*
                    Propagate match up through it hierarchy.
                    If there are Superfamily, Group, and Family models in a tree,
                    and a sequence matches F, it should inherit the S and G annotations
                */
                parents.each { parent ->
                    Signature signature = new Signature(parent, library)
                    Match promotedMatch = new Match(parent, match.evalue, match.score, match.bias, signature)
                    Location location = match.locations[0].clone()
                    location.sites.clear()
                    promotedMatch.addLocation(location)
                    matchesWithPromoted.add(promotedMatch)
                }
            }
        }

        // Merge locations from the same model together
        matches = [:]
        matchesWithPromoted.each { match ->
            Match mergedMatch = matches.get(match.modelAccession)
            if (mergedMatch) {
                mergedMatch.addLocation(match.locations[0])
            } else {
                matches[match.modelAccession] = match
            }
        }

        // Remove nested locations
        matches.each { modelAccession, match ->
            def locations = []
            match.locations
                .sort { a, b ->
                    a.start <=> b.start ?: b.end <=> a.end
                }
                .each { location ->
                    boolean isContained = locations.any { loc ->
                        loc.start <= location.start && loc.end >= location.end
                    }
                    if (!isContained) {
                        locations << location
                    }
                }

            match.locations = locations
        }

        return [seqId, matches]
    }

    def json = JsonOutput.toJson(sequences)
    new File(outputFilePath.toString()).write(json)
}

Map<String, Set<String>> getHierarchies(String filePath) {
    def hierarchies = [:].withDefault { [] as Set }
    new File(filePath).eachLine { line ->
        def nodes = line.split(/\t/).toList()
        nodes.eachWithIndex { node, idx ->
            if (idx > 0) {
                def ancestors = nodes.subList(0, idx) as Set
                hierarchies[node].addAll(ancestors)
            }
        }
    }
    return hierarchies
}

Map<String, Map<String, Match>> parseOutput(
    String outputFilePath,
    Map<String, Map> hmmerMatches
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

            assert hmmerMatches.containsKey(sequenceId)
            def seqHmmerMatches = hmmerMatches[sequenceId]
            matches[sequenceId] = parseBlock(reader, sequenceId, seqHmmerMatches)
        }
    }

    return matches
}

Map<String, Match> parseBlock(
    Reader reader,
    String sequenceId,
    Map<String, Map> seqHmmerMatches
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
            Integer hmmStart = fields[4] as Integer
            Integer hmmEnd = fields[5] as Integer
            int aliStart = fields[7] as int
            int aliEnd = fields[8] as int
            Integer envStart = fields[9] as Integer
            Integer envEnd = fields[10] as Integer
            assert seqHmmerMatches.containsKey(modelAccession)
            def hmmerMatch = seqHmmerMatches[modelAccession]
            Location loc = hmmerMatch.locations.find { it.start == aliStart 
                                                        && it.end == aliEnd 
                                                        && it.hmmStart == hmmStart 
                                                        && it.hmmEnd == hmmEnd 
                                                        && it.envelopeStart == envStart 
                                                        && it.envelopeEnd == envEnd }
            assert loc != null
            Match match = domains.computeIfAbsent(modelAccession, k -> {
                Signature signature = new Signature(modelAccession, library)
                return new Match(modelAccession, hmmerMatch.evalue, hmmerMatch.score, hmmerMatch.bias, signature)
            });
            match.addLocation(loc)
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
