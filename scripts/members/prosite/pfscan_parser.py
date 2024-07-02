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
                    return None
                if line_strip.startswith(">"):
                    seq2match = line_strip.split(":")
                    seq_id = seq2match[0].strip()
                    match_id = seq2match[1].split()[0].strip()[1:]
                    name2desc = seq2match[1].strip().split(' ', 2)
                    name = name2desc[1]
                    desc = name2desc[2]
                    if seq_id not in ips6_matches:
                        ips6_matches[seq_id] = {}
                    ips6_matches[seq_id].update({
                        match_id: {
                            "accession": match_id,
                            "name": name,
                            "description": desc,
                            "member_db": "PROSITE_PATTERNS",
                            "version": version,
                            "locations": []
                        }
                    })
                else:
                    locations = line_strip.split()
                    cigar_alignment = cigar_alignment_parser(locations[3])
                    location = {
                        "start": int(locations[0]),
                        "end": int(locations[2]),
                        "cigarAlignment": encode(cigar_alignment),
                        "alignment": locations[3]
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
