# pfam-a.hmm.path=${data.directory}/pfam/35.0/pfam_a.hmm
# pfam-a.seed.path=${data.directory}/pfam/35.0/pfam_a.seed
# pfam-clans.path=${data.directory}/pfam/35.0/pfam_clans
# pfam-a.dat.path=${data.directory}/pfam/35.0/pfam_a.dat
# hmmer3.hmmsearch.switches.pfama=-Z 17929 --cut_ga
# pfam.min.length=8


# <property name="fullPathToHmmScanBinary" value="${binary.hmmer33.hmmscan.path}"/>
# <property name="fullPathToHmmsearchBinary" value="${binary.hmmer33.hmmsearch.path}"/>
# <property name="binarySwitches" value="${hmmer3.hmmsearch.switches.pfama} ${hmmer3.hmmsearch.cpu.switch.pfama}"/>
# <property name="fullPathToHmmFile" value="${pfam-a.hmm.path}"/>



# Hmmer3SearchMatchParser


#     public void addMatch(HmmSearchRecord methodMatches, Map<String, RawProtein<PfamHmmer3RawMatch>> rawResults) throws IOException {
#         try {
#             for (SequenceMatch sequenceMatch : methodMatches.getSequenceMatches().values()) {
#                 for (DomainMatch domainMatch : sequenceMatch.getDomainMatches()) {
#                     /*
#                     Either retrieve the correct RawSequenceIdentifer, or create a new one
#                     and add it to the Map.
#                     */
#                     RawProtein<PfamHmmer3RawMatch> rawProtein = rawResults.get(sequenceMatch.getSequenceIdentifier());
#                     if (rawProtein == null) {
#                         rawProtein = new RawProtein<PfamHmmer3RawMatch>(sequenceMatch.getSequenceIdentifier());
#                         rawResults.put(sequenceMatch.getSequenceIdentifier(), rawProtein);
#                     }
#
#                     final PfamHmmer3RawMatch match = new PfamHmmer3RawMatch(
#                             sequenceMatch.getSequenceIdentifier(),
#                             methodMatches.getModelAccession(),
#                             signatureLibrary,
#                             signatureLibraryRelease,
#                             domainMatch.getAliFrom(),
#                             domainMatch.getAliTo(),
#                             sequenceMatch.getEValue(),
#                             sequenceMatch.getScore(),
#                             domainMatch.getHmmfrom(),
#                             domainMatch.getHmmto(),
#                             domainMatch.getHmmBounds(),
#                             domainMatch.getScore(),
#                             domainMatch.getEnvFrom(),
#                             domainMatch.getEnvTo(),
#                             domainMatch.getAcc(),
#                             sequenceMatch.getBias(),
#                             domainMatch.getCEvalue(),
#                             domainMatch.getIEvalue(),
#                             domainMatch.getBias()
#                     );
#                     rawProtein.addMatch(match);
#                 }



# GaValuesRetriever

# PfamHMMER3PostProcessing

# ClanFileParser