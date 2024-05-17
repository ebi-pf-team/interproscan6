import argparse
import json
import re

from pathlib import Path


SESSION_BLOCK_START_MARKER = "SESSION"
QUERY_BLOCK_START_MARKER = "QUERY"

QUERY_LINE_PATTERN = re.compile(r"^QUERY\s+(\S+)\s+(\S+)\s+(\d+)\s+(.*)$")
DOMAIN_LINE_PATTERN = re.compile(r"^(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)")
SITE_LINE_PATTERN = re.compile(r"^(\d+)\s+(\S+)\s+(\S+)\s+(.+?)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)$")
MOTIF_LINE_PATTERN = re.compile(r"^\d+\s+\w+\_\d+\s+\w+\s+\w+\s+\w+\s+\d+$")


def build_parser() -> argparse.ArgumentParser:
    """Build cmd-line argument parser"""
    parser = argparse.ArgumentParser(
        prog="cdd_match_parser",
        description="Parse the output from the CDD application, extracing matches and sites information",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "rpsblast_processed",
        type=Path,
        help="Path to rpsblast_processed file"
    )

    parser.add_argument(
        "release",
        type=str,
        help="Release version of CDD"
    )

    return parser


def parse_cdd(rpsblast_processed: Path, release: str):
    """Parse rpsblast_processed file

    :param rpsblast_processed: Path to rpsblast_processed file
    :param release: release version of CDD
    """
    matches = {}  # prot seq id: [{domain data, 'sites': {[sites]}}]
    protein_identifier = ""

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
                    sites = line.split()[-4].split(",")
                    site_info = {
                        "description": " ".join(line.split()[3:-4]),
                        "numLocations": len(sites),
                        "siteLocations": [],
                    }
                    for site in sites:
                        site_data = {
                            "start": site[1:],
                            "end": site[1:],
                            "residue": site[0]
                        }
                        try:
                            site_info["siteLocations"].append(site_data)
                        except KeyError:
                            site_info["siteLocations"] = [site_data]
                    matches[target_key][signature_accession]['locations'][-1]['sites'].append(site_info)
                    continue

    return matches


def main():
    parser = build_parser()
    args = parser.parse_args()

    matches = parse_cdd(args.rpsblast_processed, args.release)
    print(json.dumps(matches, indent=2))


if __name__ == "__main__":
    main()
