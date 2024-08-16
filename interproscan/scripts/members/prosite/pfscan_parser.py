"""
Parse the output from ps_scan.pl into the internal
IPS6 JSON structure.
"""
import json
import sys

from cigar_alignment import cigar_alignment_parser, encode


def parse(pfscan_out: str):
    with open(pfscan_out, 'r') as fh:
        ips6_matches = {}
        for line in fh:
            line_strip = line.strip()
            match_info = line_strip.split('\t')
            if line_strip:
                if len(match_info) < 9 or line.startswith("pfscanV3 is not meant to be used with a single profile"):
                    return ""
                else:
                    match_details = match_info[8].split(';')
                    level = match_details[1].strip()
                    if not level.startswith("LevelTag") or "0" not in level:  # only strong matches
                        continue
                    seq_id = match_info[0]
                    match_id = match_info[2]
                    name = match_details[0].replace('Name ', '').replace('"', '').strip()
                    alignment = match_details[2].replace('Sequence ', '').replace('"', '').replace('.', '').strip()
                    cigar_alignment = cigar_alignment_parser(alignment)
                    location = {
                        "start": int(match_info[3]),
                        "end": int(match_info[4]),
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
                            "locations": []
                        }
                    ips6_matches[seq_id][match_id]["locations"].append(location)

    return ips6_matches


def main():
    pfscan_out = sys.argv[1]  # output file from PFSCAN_RUNNER module
    output_file = sys.argv[2]  # str rep of path for internal IPS6 JSON

    ips6_matches = parse(pfscan_out)

    with open(output_file, "w") as fh:
        json.dump(ips6_matches, fh, indent=2)


if __name__ == "__main__":
    main()
