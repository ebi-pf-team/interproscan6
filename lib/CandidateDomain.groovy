class CandidateDomain {
    // Used during selection of representative domains
    Location location
    Set<Integer> residues = new HashSet<>()
    Integer representativeRank

    CandidateDomain(Location location, Integer representativeRank) {
        this.location = location
        this.representativeRank = representativeRank
    }

    Set<Integer> getResidues() {
        // Use this.@residues to access the field directly preventing a StackOverFlow from recursive calls
        this.location.fragments.each { frag ->
            this.@residues.addAll((frag.start..frag.end).toSet())
        }
        return this.@residues
    }

    List<LocationFragment> sortFragments() {
        if (this.location.fragments.size() > 1) {
            this.location.fragments.sort { a, b ->
                int comparison = a.start <=> b.start
                comparison != 0 ? comparison : a.end <=> b.end
            }
        }
        return this.location.fragments
    }
}