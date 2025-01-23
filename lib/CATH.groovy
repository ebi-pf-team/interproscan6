import HMMER3

class CATH {
    static parseAssignedFile(String filePath) {
        // For CATH-Gene3D
        def results = [:].withDefault { [] }

        new File(filePath).eachLine { line ->
            if (line[0] != "#") {
                // #domain_id cath-superfamily query-id match-id score boundaries resolved aligned-regions cond-evalue indp-evalue [comment]
                def fields = line.split("\t")
                assert fields.size() == 10 || fields.size() == 11
                String sequenceId = fields[2]
                def dom = new CathDomain(
                    fields[0],
                    fields[3],
                    "G3DSA:${fields[1]}",
                    Double.parseDouble(fields[4]),
                    Double.parseDouble(fields[9]),
                    fields[5].split(",").collect { new SimpleLocation(it) },
                    fields[6].split(",").collect { new SimpleLocation(it) },
                )

                results[sequenceId] << dom
            }        
        }

        return results
    }

    static parseResolvedFile(String filePath) {
        // For CATH-FunFam
        def results = [:].withDefault { [] }

        new File(filePath).eachLine { line ->
            if (line[0] != "#") {
                // #FIELDS query-id match-id score boundaries resolved aligned-regions cond-evalue indp-evalue
                def fields = line.split(/\s+/)
                assert fields.size() == 8
                String sequenceId = fields[0]
                String modelId = fields[1]
                String accession = "G3DSA:${modelId.replaceAll('-', ':')}"
                def dom = new CathDomain(
                    modelId,
                    modelId,
                    accession,
                    Double.parseDouble(fields[2]),
                    Double.parseDouble(fields[7]),
                    fields[3].split(",").collect { new SimpleLocation(it) },
                    fields[4].split(",").collect { new SimpleLocation(it) },
                )

                results[sequenceId] << dom
            }        
        }

        return results
    }

    static mergeWithHmmerMatches(resolvedDomains, hmmerMatches, memberDb) {
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
                    cathDomain.getResolvedStart(),
                    cathDomain.getResolvedEnd(),
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
                    domain.signature = new Signature(cathDomain.accession)
                    domain.signature.signatureLibraryRelease.library = memberDb
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
    String accession
    Double score
    Double evalue
    List<SimpleLocation> boundaries
    List<SimpleLocation> resolvedBoundaries

    CathDomain(String domainId, String matchId, String accession, 
               Double score, Double evalue, List<SimpleLocation> boundaries,
               List<SimpleLocation> resolvedBoundaries) {
        this.domainId = domainId
        this.matchId = matchId
        this.accession = accession
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

        String matchId
        if (this.matchId.startsWith("dc_")) {
            matchId = this.matchId.replace("_${this.domainId}", "")
        } else {
            matchId = this.matchId
        }

        return "${matchId}-${leftMost}-${rightMost}"
    }

    int getStart() {
        return this.boundaries*.start.min()
    }

    int getEnd() {
        return this.boundaries*.end.max()
    }

    int getResolvedStart() {
        return this.resolvedBoundaries*.start.min()
    }

    int getResolvedEnd() {
        return this.resolvedBoundaries*.end.max()
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