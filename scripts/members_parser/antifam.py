import re
import argparse
import json
from hmmer_parser import parse


def antifam_parse(hmmer_parse_result):
    # I need to improve the parse
    pass
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


def main():
    parser = argparse.ArgumentParser(
        description="Antifam parser"
    )
    parser.add_argument(
        "-seq", "--sequences", type=str, help="fasta file with sequences"
    )
    parser.add_argument("-preproc", "--preproc", type=str, help="file result of hmmer preproc")
    args = parser.parse_args()

    accession_regex = re.compile("^Accession:\\s+(ANF\\d{5})\\s*$")
    hmmer_parse_result = parse(args.preproc, accession_regex)
    #antifam_parse = ncbifam_parse(hmmer_parse_result)
    print("ANTIFAM:")
    print(json.dumps(hmmer_parse_result))


if __name__ == "__main__":
    main()
