import json
import re
import sys

from pathlib import Path


SESSION_BLOCK_START_MARKER = "SESSION"
QUERY_BLOCK_START_MARKER = "QUERY"

QUERY_LINE_PATTERN = re.compile(r"^QUERY\s+(\S+)\s+(\S+)\s+(\d+)\s+(.*)$")
DOMAIN_LINE_PATTERN = re.compile(r"^(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)")
SITE_LINE_PATTERN = re.compile(r"^(\d+)\s+(\S+)\s+(\S+)\s+(.+?)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)$")
MOTIF_LINE_PATTERN = re.compile(r"^\d+\s+\w+\_\d+\s+\w+\s+\w+\s+\w+\s+\d+$")


def parse_cdd(rpsblast_processed: Path, release: str):
    """Parse rpsblast_processed file

    :param rpsblast_processed: Path to rpsblast_processed file
    :param release: release version of CDD
    """
    matches = {}  # prot seq id: [{domain data, 'sites': {[sites]}}]
    protein_identifier = ""
    signature_accession = ""

    with open(rpsblast_processed, "r") as fh:
        for line in fh.readlines():

            if line.startswith(SESSION_BLOCK_START_MARKER):
                # #SESSION	<session-ordinal>	<program>	<database>	<score-matrix>	<evalue-threshold>
                # e.g. SESSION	1	blastp	2.15.0+	.../InterProScan6/data/cdd/3.20/db/Cdd_NCBI	BLOSUM62	0.01
                protein_identifier = line.split()[1].strip()

            elif line.startswith(QUERY_BLOCK_START_MARKER):
                # #QUERY	<query-id>	<seq-type>	<seq-length>	<definition-line>
                # e.g. QUERY BLOCK: QUERY      Query_32        Peptide 576     A0A095XSD9_9FIRM 1-576
                query_line = QUERY_LINE_PATTERN.match(line)
                if query_line:
                    definition_line = query_line.group(4).strip()
                    target_key = definition_line.split("|")[1] if definition_line.startswith("sp|") else definition_line.split()[0]
                    # target_key = protein seq id from input file

                    if target_key not in matches:
                        matches[target_key] = {}

            elif line.startswith(protein_identifier):
                _line = DOMAIN_LINE_PATTERN.match(line)
                if _line:
                    # #<session-ordinal>	<query-id[readingframe]>	<hit-type>	<PSSM-ID>	<from>	<to>	<E-Value>	<bitscore>	<accession>	<short-name>	<incomplete>	<superfamily PSSM-ID>
                    # e.g. 1	Query_1	Non-specific	238121	235	438	1.87425e-09	56.9596	cd00200	WD40	N	453027
                    if line.split()[2].upper() == "SPECIFIC":  # filter for only matches of hit-type SPECIFC
                        signature_accession = _line.group(9)

                        domain_info = {
                            "start": int(_line.group(5)),
                            "end": int(_line.group(6)),
                            "representative": "",
                            "evalue": float(_line.group(7)),
                            "score": float(_line.group(8)),
                            "postProcessed": "",
                            "sites": [],
                        }

                        signature_info = {
                            "accession": signature_accession,
                            "name": _line.group(10),
                            "member_db": "CDD",
                            "version": release,
                            "model-ac": signature_accession,
                        }

                        if signature_accession not in matches[target_key]:
                            matches[target_key][signature_accession] = signature_info
                            matches[target_key][signature_accession]['locations'] = [domain_info]
                        else:
                            matches[target_key][signature_accession]["locations"].append(domain_info)

                    continue

                _line = SITE_LINE_PATTERN.match(line)
                if _line:
                    if signature_accession in matches[target_key]:
                        sites = line.split()[-4].split(",")

                        # find the domain that contains these sites
                        positions = [int(_[1:].split("-")[0]) for _ in sites]
                        earliest_start, latest_end = min(positions), max(positions)
                        domain_index = None
                        for i in range(len(matches[target_key][signature_accession]['locations'])):
                            if (
                                matches[target_key][signature_accession]['locations'][i]["start"] <= earliest_start
                            ) and (
                                matches[target_key][signature_accession]['locations'][i]["end"] >= latest_end
                            ):
                                domain_index = i
                                break

                        if domain_index is not None:  # domain index could by 0, so don't use `if domain_index:`
                            site_info = {
                                "description": " ".join(line.split()[3:-4]),
                                "numLocations": len(sites),
                                "siteLocations": [],
                            }
                            site_info["siteLocations"] = []
                            for site in sites:
                                site_info["siteLocations"].append(
                                    {
                                        "start": int(site[1:]),
                                        "end": int(site[1:]),
                                        "residue": site[0]
                                    }
                                )
                            matches[target_key][signature_accession]['locations'][-1]['sites'].append(site_info)
                    continue

    return matches


def main():
    """CL input:
    0. Str repr of path to the rpsblast processed file
    1. release
    2. Str repr of path for the output file
    """
    args = sys.argv[1:]
    matches = parse_cdd(args[0], args[1])
    with open(args[2], "w") as fh:
        json.dump(matches, fh)


if __name__ == "__main__":
    main()