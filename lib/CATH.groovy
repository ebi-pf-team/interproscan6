import HMMER3

class CATH {
    static parseAssignedFile(String filePath) {
        def results = [:].withDefault { [] }

        new File(filePath).eachLine { line ->
            if (line[0] != "#") {
                def fields = line.split("\t")
                assert fields.size() == 10
                String sequenceId = fields[2]
                def dom = new CathDomain(
                    fields[0],                                                  // Domain ID
                    fields[3],                                                  // Match ID
                    fields[1],                                                  // CATH ID
                    Double.parseDouble(fields[4]),                              // Score
                    Double.parseDouble(fields[9]),                              // i-evalue
                    fields[5].split(",").collect { new SimpleLocation(it) },    // Boundaries
                    fields[6].split(",").collect { new SimpleLocation(it) },    // Resolved
                )

                results[sequenceId] << dom
            }        
        }

        return results
    }

    static mergeWithHmmerMatches(resolvedDomains, hmmerMatches) {
        def hmmerDomains = HMMER3.splitByLocation(hmmerMatches)
        def results = [:].withDefault { [:] }
        resolvedDomains.each { sequenceId, domains ->
            def thiSeqHmmerDomains = hmmerDomains.get(sequenceId)
            assert thiSeqHmmerDomains != null

            domains.each { cathDomain ->
                def key = cathDomain.getKey()
                def hmmerDomain = thiSeqHmmerDomains.get(key)
                assert hmmerDomain != null

                def boundaries = cathDomain.resolvedBoundaries
                def fragments = []
                if (boundaries.size() > 1) {
                    cathDomain.resolvedBoundaries.eachWithIndex { x, i ->
                        LocationFragment fragment
                        if (i == 0) {
                            fragment = new LocationFragment(x.start, x.end, "C_TERMINAL_DISC")
                        } else if (i + 1 < boundaries.size()) {
                            fragment = new LocationFragment(x.start, x.end, "NC_TERMINAL_DISC")
                        } else {
                            fragment = new LocationFragment(x.start, x.end, "N_TERMINAL_DISC")
                        }
                        fragments.add(fragment)
                    }
                } else {
                    LocationFragment fragment = new LocationFragment(boundaries[0].start, boundaries[0].end, "CONTINUOUS")
                    fragments.add(fragment)
                }

                Location location = new Location(
                    cathDomain.getStart(),
                    cathDomain.getEnd(),
                    hmmerDomain.locations[0].hmmStart,
                    hmmerDomain.locations[0].hmmEnd,
                    hmmerDomain.locations[0].hmmLength,
                    hmmerDomain.locations[0].hmmBounds,
                    hmmerDomain.locations[0].envelopeStart,
                    hmmerDomain.locations[0].envelopeEnd,
                    cathDomain.evalue,
                    cathDomain.score,
                    null,
                    fragments
                )

                def sequenceDomains = results[sequenceId]
                def domId = cathDomain.domainId
                if (sequenceDomains.containsKey(domId)) {
                    sequenceDomains[domId].addLocation(location)
                } else {
                    Match domain = new Match(
                        domId, 
                        hmmerDomain.evalue,
                        hmmerDomain.score, 
                        hmmerDomain.bias
                    )
                    domain.signature = new Signature(cathDomain.supfamId)
                    domain.addLocation(location)
                    sequenceDomains[domId] = domain
                }
            }
        }
        return results
    }
}

class CathDomain {
    String domainId
    String matchId
    String supfamId
    Double score
    Double evalue
    List<SimpleLocation> boundaries
    List<SimpleLocation> resolvedBoundaries

    CathDomain(String domainId, String matchId, String supfamId, 
               Double score, Double evalue, List<SimpleLocation> boundaries,
               List<SimpleLocation> resolvedBoundaries) {
        this.domainId = domainId
        this.matchId = matchId
        this.supfamId = "G3DSA:${supfamId}"
        this.score = score
        this.evalue = evalue       
        this.boundaries = this.sortLocations(boundaries)
        this.resolvedBoundaries = this.sortLocations(resolvedBoundaries)
    }

    static List<SimpleLocation> sortLocations(List<SimpleLocation> locations) {
        return locations.sort { a, b ->
            a.start <=> b.start ?: a.end <=> b.end
        }
    }

    String getKey() {
        int leftMost = this.getStart()
        int rightMost = this.getEnd()
        return "${this.matchId}-${leftMost}-${rightMost}"
    }

    int getStart() {
        return this.boundaries*.start.min()
    }

    int getEnd() {
        return this.boundaries*.end.max()
    }
}

class SimpleLocation {
    int start
    int end

    SimpleLocation(int start, int end) {
        this.start = start
        this.end = end
    }

    SimpleLocation(String range) {
        String[] fields = range.split("-")
        assert fields.size() == 2
        this.start = fields[0].toInteger()
        this.end = fields[1].toInteger()
    }
}