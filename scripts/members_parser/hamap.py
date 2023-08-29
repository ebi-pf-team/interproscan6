import re

project_dir = "/Users/lcf/interproscan6"

pfsearch_wrapperpath = f"{project_dir}/bin/prosite/pfsearch_wrapper.py"
profile_models_path = f"{project_dir}/data/hamap/2023_01/hamap.prf"
profile_models_dir = f"{project_dir}/data/hamap/2023_01/profiles"
hmm_path = f"{project_dir}/data/hamap/2023_01/hamap.hmm.lib"
# hmmer3_hmmsearch_switches.hmmfilter=-E 100 --domE 100 --incE 100 --incdomE 100
# pfsearchv3_hamap_binary.switches=-f -o 7


def check_parser():
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


# <property name="fullPathToHmmScanBinary" value="${binary.hmmer33.hmmscan.path}"/>
# <property name="binarySwitches" value="${hmmer3.hmmsearch.switches.hmmfilter} ${hmmer3.hmmsearch.cpu.switch.hmmfilter}"/>
# <property name="fullPathToHmmFile" value="${hamap.hmm.path}"/>
# <property name="useTbloutFormat" value="true"/>
# <property name="outputFileNameTbloutTemplate" ref="rawAnalysisOutputTbloutFileTemplate"/>


if __name__ == "__main__":
    pass
