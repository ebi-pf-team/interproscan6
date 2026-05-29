process SEARCH_SFLD {
    label     'mem_min'
    label     'time_veryshort'
    label     'dynamic'
    container 'interpro/sfld:ebiftp'

    input:
    tuple val(meta), path(fasta)
    tuple path(hmm), path(sites)

    output:
    tuple val(meta), path("hmmsearch.out"), path("sfld.tsv")

    script:
    """
    hmmsearch \
        -Z 378 --acc \
        --cut_ga \
        --cpu ${task.cpus} \
        -o hmmsearch.out \
        --domtblout hmmsearch.tab \
        -A hmmsearch.sto \
        ${hmm} ${fasta}

    sfld_postprocess \
        --alignment hmmsearch.sto \
        --dom hmmsearch.tab \
        --hmmer-out hmmsearch.out \
        --site-info ${sites} \
        --output sfld.tsv
    """
}

process PARSE_SFLD {
    label    'mem_low'
    label    'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(hmmsearch_out), val(postprocess_out)
    val hierarchydb

    output:
    tuple val(meta), path("sfld.json")

    exec:
    def hmmerMatches = uk.ac.ebi.interpro.HMMER3.parseOutput(hmmsearch_out, "SFLD")
    def sequences = uk.ac.ebi.interpro.SFLD.parseOutput(postprocess_out, hmmerMatches)
    def hierarchies = uk.ac.ebi.interpro.SFLD.getHierarchies(hierarchydb)
    def library = new uk.ac.ebi.interpro.SignatureLibraryRelease("SFLD", null)

    sequences = sequences.collectEntries { seqId, matches -> 
        // Flatten matches (one location per match)
        matches = matches.collectMany { _key, match ->
            return match.locations.collect { location ->
                def signature = new uk.ac.ebi.interpro.Signature(match.modelAccession, library)
                def newMatch = new uk.ac.ebi.interpro.Match(match.modelAccession, match.evalue, match.score, match.bias, signature)
                newMatch.addLocation(location.clone())
                return newMatch
            }
        }

        if (matches.size() > 1) {
            // Remove overlapping matches from models in the same hierarchy (keep the most specific)
            def ignored = [] as Set
            matches = matches
                .collect { match ->
                    if (match.modelAccession.startsWith("SFLDF")) {
                        // Highly specific model (family): always keep
                        return match
                    }

                    def matchAccession = match.modelAccession
                    def l1 = match.locations[0]
                    def ignore = matches.any { otherMatch ->
                        if (matchAccession == otherMatch.modelAccession) {
                            // Two matches/locations from the same mode: keep
                            return false
                        }
                        if (ignored.contains(otherMatch)) {
                            // otherMatch has previously been flagged as to be ignored
                            return false
                        }

                        def l2 = otherMatch.locations[0]
                        if (l1.start > l2.end || l2.start > l1.end) {
                            return false
                        }

                        /*
                            The current match overlaps a match from a more specific model:
                            we want to keep the most specific match, so we ignore the current match
                        */
                        def otherParents = hierarchies.get(otherMatch.modelAccession)
                        return otherParents != null && otherParents.contains(matchAccession)
                    }

                    if (ignore) {
                        // Mark this match as to be ignored
                        ignored.add(match)
                        return null
                    }

                    return match
                }
                .findAll { m -> m != null }
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
                    def signature = new uk.ac.ebi.interpro.Signature(parent, library)
                    def promotedMatch = new uk.ac.ebi.interpro.Match(parent, match.evalue, match.score, match.bias, signature)
                    def location = match.locations[0].clone()
                    location.sites.clear()
                    promotedMatch.addLocation(location)
                    matchesWithPromoted.add(promotedMatch)
                }
            }
        }

        // Merge locations from the same model together
        matches = [:]
        matchesWithPromoted.each { match ->
            def mergedMatch = matches.get(match.modelAccession)
            if (mergedMatch) {
                mergedMatch.addLocation(match.locations[0])
            } else {
                matches[match.modelAccession] = match
            }
        }

        // Remove nested locations
        matches.each { _modelAccession, match ->
            def locations = []
            match.locations
                .sort { a, b ->
                    a.start <=> b.start ?: b.end <=> a.end
                }
                .each { location ->
                    def isContained = locations.any { loc ->
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

    def filepath = task.workDir.resolve("sfld.json")
    filepath.text = groovy.json.JsonOutput.toJson(sequences)
}
