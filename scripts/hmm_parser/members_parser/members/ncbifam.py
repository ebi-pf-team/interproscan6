import re


def ncbifam_parse(hmmer_parse_result):
    # I need to improve the parse
    pass
    # protected NCBIfamRawMatch createMatch(final String signatureLibraryRelease,
    #                                       final HmmSearchRecord hmmSearchRecord,
    #                                       final SequenceMatch sequenceMatch,
    #                                       final DomainMatch domainMatch) {
    #     return new NCBIfamRawMatch(
    #             sequenceMatch.getSequenceIdentifier(),
    #             hmmSearchRecord.getModelAccession(),
    #             SignatureLibrary.NCBIFAM,
    #             signatureLibraryRelease,
    #             domainMatch.getAliFrom(),
    #             domainMatch.getAliTo(),
    #             sequenceMatch.getEValue(),
    #             sequenceMatch.getScore(),
    #             domainMatch.getHmmfrom(),
    #             domainMatch.getHmmto(),
    #             domainMatch.getHmmBounds(),
    #             domainMatch.getScore(),
    #             domainMatch.getEnvFrom(),
    #             domainMatch.getEnvTo(),
    #             domainMatch.getAcc(),
    #             sequenceMatch.getBias(),
    #             domainMatch.getCEvalue(),
    #             domainMatch.getIEvalue(),
    #             domainMatch.getBias()
    #     );
    # }


