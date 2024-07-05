"""
Parse the output from ps_scan.pl into the internal
IPS6 JSON structure.
"""
import json
import sys

from cigar_alignment import cigar_alignment_parser, encode


def parse(pfscan_out: str, version: str):
    with open(pfscan_out, 'r') as fh:
        ips6_matches = {}
        for line in fh:
            line_strip = line.strip()
            if line_strip:
                if line.startswith("pfscanV3 is not meant to be used with a single profile"):
                    return ""
                parts = line.split()
                if len(line.split()) > 9:
                    seq_id = parts[0]
                    match_id = parts[2]
                    name = parts[9].strip().replace('SequenceDescription "', '').replace('"', '')
                    alignment = parts[15].replace('"', '').replace('.', '')
                    cigar_alignment = cigar_alignment_parser(alignment)
                    location = {
                        "start": int(parts[3]),
                        "end": int(parts[4]),
                        "representative": "false",
                        "level": "STRONG",
                        "cigarAlignment": encode(cigar_alignment),
                        "alignment": alignment
                    }
                    if seq_id not in ips6_matches:
                        ips6_matches[seq_id] = {}
                    if match_id not in ips6_matches[seq_id]:
                        ips6_matches[seq_id][match_id] = {
                            "accession": match_id,
                            "name": name,
                            "description": "-",
                            "member_db": "PROSITE_PATTERNS",
                            "version": version,
                            "locations": []
                        }
                    ips6_matches[seq_id][match_id]["locations"].append(location)

    return ips6_matches


def main():
    pfscan_out = sys.argv[1]  # output file from PFSCAN_RUNNER module
    output_file = sys.argv[2]  # str rep of path for internal IPS6 JSON

    version = pfscan_out.split("._.")[0]

    ips6_matches = parse(pfscan_out, version)

    with open(output_file, "w") as fh:
        json.dump(ips6_matches, fh, indent=2)


if __name__ == "__main__":
    main()
