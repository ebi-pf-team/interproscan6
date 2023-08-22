import re


def parser():
    Accession = re.compile("^Accession:\\s+(ANF\\d{5})\\s*$")
    # protected AntiFamHmmer3RawMatch createMatch(final String signatureLibraryRelease,
    #                                             final HmmSearchRecord hmmSearchRecord,
    #                                             final SequenceMatch sequenceMatch,
    #                                             final DomainMatch domainMatch) {
    #     return new AntiFamHmmer3RawMatch(
    #             sequenceMatch.getSequenceIdentifier(),
    #             hmmSearchRecord.getModelAccession(),
    #             SignatureLibrary.ANTIFAM,
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


if __name__ == '__main__':
    parser()
