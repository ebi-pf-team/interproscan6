import stockholm_parser


# Post-processes raw results for Pfam HMMER3 in the batch requested.
#      *
#      * @param proteinIdToRawMatchMap being a Map of protein IDs to a List of raw matches
#      * @return a Map of proteinIds to a List of filtered matches.

def process_matches(raw_protein_unfiltered, seed_matches, filtered_matches, clan_data):
    # evalue ASC e score DESC
    unfiltered_by_evalue = sorted(
        raw_protein_unfiltered,
        key=lambda m: (m['evalue'], -m['score'])
    )
    for candidate_match in unfiltered_by_evalue:
        if candidate_match not in seed_matches:
            candidate_match_clan = clan_data.get(candidate_match['model_id'])
            passes = True
            if candidate_match_clan:
                for match in filtered_matches:
                    passed_match_clan = clan_data.get(match['model_id'])
                    # same cla?
                    if candidate_match_clan == passed_match_clan:
                        # overlap?
                        if matches_overlap(candidate_match, match):
                            if not matches_are_nested(candidate_match, match):
                                passes = False
                                break
            if passes:
                filtered_matches.append(candidate_match)


def matches_overlap(one, two):
    return max(one['location_start'], two['location_start']) <= min(one['location_end'], two['location_end'])


def regions_overlap(start_region_one, end_region_one, start_region_two, end_region_two):
    return max(start_region_one, start_region_two) <= min(end_region_one, end_region_two)


def matches_are_nested(one, two, clan_data):
    one_model = clan_data.get_model_by_model_accession(one['model_id'])
    two_model = clan_data.get_model_by_model_accession(two['model_id'])
    if one_model is None or two_model is None:
        return False
    return one_model.is_nested_in(two_model) or two_model.is_nested_in(one_model)


def main():
    # args = sys.argv[1:]
    #
    # clan_data_parsed = args[0]



if __name__ == '__main__':
    main()





    # private PfamClanData clanData; # o raw do clan?
    #
    # private ClanFileParser clanFileParser; # parser do clan
    #
    # private SeedAlignmentDataRetriever seedAlignmentDataRetriever; # matches do run hmmer, acho

    # /**
    #  * @param seedAlignmentDataRetriever to retrieve seed alignment data for
    #  *                                   a range of proteins.
    #  */
    #
    # /**
    #  * Post-processes raw results for Pfam HMMER3 in the batch requested.
    #  *
    #  * @param proteinIdToRawMatchMap being a Map of protein IDs to a List of raw matches
    #  * @return a Map of proteinIds to a List of filtered matches.
    #  */


    # Also includes additional code to ensure seed alignments are included as matches, regardless of
#      * score.
#      *
#      * @param rawProteinUnfiltered being a List of the raw matches to filter
#      * @param seedAlignments       being a Collection of SeedAlignment objects, to check for matches to
#      *                             methods where this protein was part of the seed alignment.
#      * @return a List of filtered matches.


    # I need to understand better this part
    private RawProtein processProtein(final RawProtein<PfamHmmer3RawMatch> rawProteinUnfiltered, final Map<String, Set<String>> nestedModelsMap, final List<SeedAlignment> seedAlignments) {
        RawProtein<PfamHmmer3RawMatch> filteredMatches = new RawProtein<PfamHmmer3RawMatch>(rawProteinUnfiltered.getProteinIdentifier());
        RawProtein<PfamHmmer3RawMatch> filteredRawProtein = new RawProtein<PfamHmmer3RawMatch>(rawProteinUnfiltered.getProteinIdentifier());

        // First of all, place any rawProteinUnfiltered to methods for which this protein was a seed
        // into the filteredMatches collection.
        final Set<PfamHmmer3RawMatch> seedMatches = new HashSet<PfamHmmer3RawMatch>();

        # remove??
        # // seedAlignments is ALWAYS null
        # if (seedAlignments != null) {        // TODO This check can be removed, once the seed alignment stuff has been sorted.
        #     Utilities.verboseLog(localVerboseLevel,"seedAlignments count:" + seedAlignments.size());
        #     for (final SeedAlignment seedAlignment : seedAlignments) {
        #         for (final PfamHmmer3RawMatch candidateMatch : rawProteinUnfiltered.getMatches()) {
        #
        #             if (!seedMatches.contains(candidateMatch)) {
        #                 if (seedAlignment.getModelAccession().equals(candidateMatch.getModelId()) &&
        #                         seedAlignment.getAlignmentStart() <= candidateMatch.getLocationStart() &&
        #                         seedAlignment.getAlignmentEnd() >= candidateMatch.getLocationEnd()) {
        #                     // Found a match to a seed, where the coordinates fall within the seed alignment.
        #                     // Add it directly to the filtered rawProteinUnfiltered...
        #                     Utilities.verboseLog(localVerboseLevel,"found match to a seed - candidateMatch and seedMatch: " + candidateMatch);
        #                     filteredMatches.addMatch(candidateMatch);
        #                     seedMatches.add(candidateMatch);
        #                 }
        #             }
        #         }
        #     }
        # }

        // Then iterate over the non-seed raw rawProteinUnfiltered, sorted in order ievalue ASC score DESC
        final Set<PfamHmmer3RawMatch> unfilteredByEvalue = new TreeSet<PfamHmmer3RawMatch>(rawProteinUnfiltered.getMatches());

        for (final RawMatch rawMatch : unfilteredByEvalue) {
            final PfamHmmer3RawMatch candidateMatch = (PfamHmmer3RawMatch) rawMatch;
            Utilities.verboseLog(localVerboseLevel,"consider match - candidateMatch: " + candidateMatch);
            if (!seedMatches.contains(candidateMatch)) {
                final PfamClan candidateMatchClan = clanData.getClanByModelAccession(candidateMatch.getModelId());

                boolean passes = true;   // Optimistic algorithm!
                Utilities.verboseLog(localVerboseLevel,"candidateMatchClan: " + candidateMatchClan);
                if (candidateMatchClan != null) {
                    // Iterate over the filtered rawProteinUnfiltered (so far) to check for passes
                    for (final PfamHmmer3RawMatch match : filteredMatches.getMatches()) {
                        final PfamClan passedMatchClan = clanData.getClanByModelAccession(match.getModelId());
                        // Are both the candidate and the passedMatch in the same clan?
                        if (candidateMatchClan.equals(passedMatchClan)) {
                            // Both in the same clan, so check for overlap.  If they overlap
                            // and are NOT nested, then set passes to false and break out of the inner for loop.
                            if (matchesOverlap(candidateMatch, match)) {
                                if (!matchesAreNested(candidateMatch, match)) {
                                    passes = false;
                                    break;  // out of loop over filtered rawProteinUnfiltered.
                                } else {
                                    Utilities.verboseLog(localVerboseLevel,"nested match: candidateMatch - " + candidateMatch
                                            + " other match:- " + match);
                                }
                            }
                        }
                    }
                }

                if (passes) {
                    // Add filtered match to collection
                    filteredMatches.addMatch(candidateMatch);
                }
            }
        }

        for (PfamHmmer3RawMatch pfamHmmer3RawMatch : filteredMatches.getMatches()) {
            String modelId = pfamHmmer3RawMatch.getModelId();

            Set<String> nestedModels = nestedModelsMap.get(modelId);
            Utilities.verboseLog(localVerboseLevel,"nestedModels: " + nestedModels);
            if (nestedModels != null && ! nestedModels.isEmpty()) {
                final UUID splitGroup = UUID.randomUUID();
                pfamHmmer3RawMatch.setSplitGroup(splitGroup);
                //get new regions
                List<Hmmer3Match.Hmmer3Location.Hmmer3LocationFragment> locationFragments = new ArrayList<>();
                int nestedFragments = 0;
                for (PfamHmmer3RawMatch rawMatch : filteredMatches.getMatches()) {
                    if (nestedModels.contains(rawMatch.getModelId()) &&
                            (matchesOverlap(rawMatch, pfamHmmer3RawMatch))) {
                        locationFragments.add(new Hmmer3Match.Hmmer3Location.Hmmer3LocationFragment(
                                rawMatch.getLocationStart(), rawMatch.getLocationEnd()));
                        nestedFragments ++;
                    }
                }
                Utilities.verboseLog(localVerboseLevel,"locationFragments to consider:  (# " + nestedFragments + ")" + locationFragments.toString());
                //the following is for testing only should be removed in the main code later
//                locationFragments.add(new Hmmer3Match.Hmmer3Location.Hmmer3LocationFragment(
//                        380, 395));
                //sort these according to the start and stop positions
                Collections.sort(locationFragments);

                DCStatus fragmentDCStatus = DCStatus.CONTINUOUS;

                List<PfamHmmer3RawMatch> rawDiscontinuousMatches  = new ArrayList<>();
                rawDiscontinuousMatches.add(pfamHmmer3RawMatch);
                if (nestedFragments > 1){
                    Utilities.verboseLog(localVerboseLevel,"nestedFragments > 1 require special investigation ");
                }
                for (Hmmer3Match.Hmmer3Location.Hmmer3LocationFragment fragment : locationFragments) {
                    List<PfamHmmer3RawMatch> newMatchesFromFragment  = new ArrayList<>();
                    for (PfamHmmer3RawMatch rawDiscontinuousMatch: rawDiscontinuousMatches) {
                        Utilities.verboseLog(localVerboseLevel,"rawDiscontinuousMatch to consider: " + rawDiscontinuousMatch.toString());
                        int newLocationStart = rawDiscontinuousMatch.getLocationStart();
                        int newLocationEnd = rawDiscontinuousMatch.getLocationEnd();
                        int finalLocationEnd = rawDiscontinuousMatch.getLocationEnd();
                        if (! regionsOverlap(newLocationStart,newLocationEnd, fragment.getStart(), fragment.getEnd())){
                            newMatchesFromFragment.add(rawDiscontinuousMatch);  // we add this match as previously processed
                            continue;
                        }
                        if (fragment.getStart() <= newLocationStart && fragment.getEnd() >= newLocationEnd){
                            fragmentDCStatus = DCStatus.NC_TERMINAL_DISC;
                            rawDiscontinuousMatch.setLocFragmentDCStatus(fragmentDCStatus.getSymbol());
                            newMatchesFromFragment.add(rawDiscontinuousMatch);
                            continue;
                        }

                        if(fragmentDCStatus ==  DCStatus.CONTINUOUS){
                            fragmentDCStatus = null;
                        }
                        boolean twoAtualRegions = false;
                        Utilities.verboseLog(localVerboseLevel,"region to consider: " + fragment.toString());
                        if (fragment.getStart() <= newLocationStart) {
                            newLocationStart = fragment.getEnd() + 1;
                            fragmentDCStatus = DCStatus.N_TERMINAL_DISC;
                        } else if (fragment.getEnd() >= newLocationEnd) {
                            newLocationEnd = fragment.getStart() - 1;
                            fragmentDCStatus = DCStatus.getNewDCStatus(fragmentDCStatus, DCStatus.C_TERMINAL_DISC);
                        } else if (fragment.getStart() > newLocationStart && fragment.getEnd() < newLocationEnd) {
                            //we have two new fragments
                            newLocationEnd = fragment.getStart() - 1;
                            twoAtualRegions = true;
                            fragmentDCStatus = DCStatus.getNewDCStatus(fragmentDCStatus,  DCStatus.C_TERMINAL_DISC);
                        }
                        Utilities.verboseLog(localVerboseLevel,"New Region: " + newLocationStart + "-" + newLocationEnd);
                        PfamHmmer3RawMatch pfMatchRegionOne = getTempPfamHmmer3RawMatch(pfamHmmer3RawMatch, newLocationStart, newLocationEnd, fragmentDCStatus);
                        pfMatchRegionOne.setSplitGroup(splitGroup);
                        pfMatchRegionOne.setLocFragmentDCStatus(fragmentDCStatus.getSymbol());
                        newMatchesFromFragment.add(pfMatchRegionOne);
                        newLocationStart = fragment.getEnd() + 1;
                        Utilities.verboseLog(localVerboseLevel," New Match for Region One  :" + pfMatchRegionOne.toString());
                        if (twoAtualRegions) {
                            //deal with final region
                            fragmentDCStatus = DCStatus.N_TERMINAL_DISC;
                            Utilities.verboseLog(localVerboseLevel,"The Last new Region: " + newLocationStart + "-" + finalLocationEnd);
                            PfamHmmer3RawMatch pfMatchRegionTwo = getTempPfamHmmer3RawMatch(pfamHmmer3RawMatch, newLocationStart, finalLocationEnd, fragmentDCStatus);
                            pfMatchRegionTwo.setSplitGroup(splitGroup);
                            pfMatchRegionTwo.setLocFragmentDCStatus(fragmentDCStatus.getSymbol());
                            newMatchesFromFragment.add(pfMatchRegionTwo);
                            Utilities.verboseLog(localVerboseLevel," New Match for Region Two :" + pfMatchRegionTwo.toString());
                        }
                    }
                    rawDiscontinuousMatches = newMatchesFromFragment;
                }
                //now add the processed discontinuous matches for further post processing or filtering into actual matches
                for (PfamHmmer3RawMatch rawDiscontinuousMatch: rawDiscontinuousMatches) {
                    int matchLength = rawDiscontinuousMatch.getLocationEnd() - rawDiscontinuousMatch.getLocationStart() + 1;
                    if (matchLength >= this.getMinMatchLength()) {
                        filteredRawProtein.addMatch(rawDiscontinuousMatch);
                    }
                }
            } else if ((pfamHmmer3RawMatch.getLocationEnd() - pfamHmmer3RawMatch.getLocationStart() + 1) >= this.getMinMatchLength()) {
                filteredRawProtein.addMatch(pfamHmmer3RawMatch);
            }
        }

        return filteredRawProtein;
    }


    # Just matches parsed (I need to check if I have all info)
    # private PfamHmmer3RawMatch getTempPfamHmmer3RawMatch(PfamHmmer3RawMatch rawMatch, int start, int end, DCStatus dcStatus) {
    #     final PfamHmmer3RawMatch match = new PfamHmmer3RawMatch(
    #             rawMatch.getSequenceIdentifier(),
    #             rawMatch.getModelId(),
    #             rawMatch.getSignatureLibrary(),
    #             rawMatch.getSignatureLibraryRelease(),
    #             start,
    #             end,
    #             rawMatch.getEvalue(),
    #             rawMatch.getScore(),
    #             rawMatch.getHmmStart(),
    #             rawMatch.getHmmEnd(),
    #             rawMatch.getHmmBounds(),
    #             rawMatch.getScore(),
    #             rawMatch.getEnvelopeStart(),
    #             rawMatch.getEnvelopeEnd(),
    #             rawMatch.getExpectedAccuracy(),
    #             rawMatch.getFullSequenceBias(),
    #             rawMatch.getDomainCeValue(),
    #             rawMatch.getDomainIeValue(),
    #             rawMatch.getDomainBias()
    #     );
    #     match.setLocFragmentDCStatus(dcStatus.getSymbol());
    #
    #     return match;
    # }
