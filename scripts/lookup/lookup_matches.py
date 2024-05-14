import json
import sys
import xml.etree.ElementTree as ET

import requests


def match_lookup(matches_checked: list, url: str) -> str:
    url_input = ','.join(matches_checked)
    matches = requests.get(f"{url}?md5={url_input}")
    return matches.text


def parse_match(match_data: str, applications: list, md52seq_id: dict) -> dict:
    tree = ET.fromstring(match_data)
    matches = {}
    hmm_bound_pattern = {"[]": "Complete", "[.": "N-terminal complete", ".]": "C-terminal complete", "..": "Incomplete"}

    for match in tree.findall(".//match"):
        for hit in match.findall("hit"):
            hit_data = hit.text.split(',')
            hit_appl = hit_data[0]
            if hit_appl in applications:
                protein_md5 = match.find("proteinMD5").text
                if protein_md5.lower() in md52seq_id:
                    target_key = md52seq_id[protein_md5.lower()]
                else:
                    target_key = protein_md5

                accession = hit_data[2]
                post_processed = False
                if hit_appl == "gene3d" or hit_appl == "pfam":
                    post_processed = True

                signature = {
                    "accession": accession,
                    "name": "",
                    "description": "",
                    "evalue": float(hit_data[16]),
                    "score": float(hit_data[7]),
                    "bias": float(hit_data[8]),
                    "version": hit_data[1],
                    "member_db": hit_appl
                }

                location = {
                    "start": int(hit_data[4]),
                    "end": int(hit_data[5]),
                    "representative": "",
                    "hmmStart": int(hit_data[10]),
                    "hmmEnd": int(hit_data[11]),
                    "hmmLength": int(hit_data[12]),
                    "hmmBounds": hit_data[9],
                    "hmmBoundsParsed": hmm_bound_pattern[hit_data[9]],
                    "evalue": float(hit_data[16]),
                    "score": hit_data[15],
                    "envelopeStart": int(hit_data[13]),
                    "envelopeEnd": int(hit_data[14]),
                    "postProcessed": post_processed,
                    "cigar_alignment": hit_data[6]
                }

                if target_key not in matches:
                    matches[target_key] = {}

                if accession not in matches[target_key]:
                    matches[target_key][accession] = signature
                    matches[target_key][accession]["locations"] = [location]
                else:
                    matches[target_key][accession]["locations"].append(location)

    return matches


def main():
    args = sys.argv[1:]

    applications = list(map(lambda x: x.upper(), args[1].split(',')))

    with open(args[0], 'r') as md5_data:
        checked_data = json.load(md5_data)
    matches = checked_data["matches"]
    seq_info = checked_data["sequences_info"]

    md52seq_id = {}
    for seq_id, match in seq_info.items():
        md52seq_id[match[-2]] = seq_id

    match_results = match_lookup(matches, args[2])
    match_parsed = parse_match(match_results, applications, md52seq_id)

    print(json.dumps(match_parsed, indent=2))


if __name__ == "__main__":
    main()
