class CandidateLocation {
    // Used during selection of representative locations
    Location location
    Set<Integer> residues = new HashSet<>()
    Integer representativeRank

    CandidateLocation(Location location, Integer representativeRank) {
        this.location = location
        this.representativeRank = representativeRank
        populateResidues()
        sortFragments()
    }

    private void populateResidues() { // Note: the private kywd is broken in groovy, so it's more an annotation
        // Use this.@residues to access the field directly preventing a StackOverFlow from recursive calls
        this.location.fragments.each { frag ->
            this.@residues.addAll((frag.start..frag.end).toSet())
        }
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