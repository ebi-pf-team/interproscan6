import json
import sys
import xml.etree.ElementTree as ET

import urllib.request


def match_lookup(matches_checked: list, url: str) -> str:
    url_input = ','.join(matches_checked)
    matches = urllib.request.urlopen(f"{url}?md5={url_input}")
    return matches.read().decode('utf-8')


def parse_match(match_data: str, applications: list, md52seq_id: dict) -> dict:
    tree = ET.fromstring(match_data)
    matches = {}
    hmm_bound_pattern = {"[]": "COMPLETE", "[.": "N_TERMINAL_COMPLETE", ".]": "C_TERMINAL_COMPLETE", "..": "INCOMPLETE"}

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
                post_processed = "false"
                if hit_appl == "gene3d" or hit_appl == "pfam":
                    post_processed = "true"
                try:
                    hmm_bounds = hmm_bound_pattern[hit_data[9]]
                except KeyError:
                    hmm_bounds = ""

                signature = {
                    "accession": accession,
                    "model-ac": hit_data[3],
                    "name": "",
                    "description": "",
                    "evalue": float(hit_data[8]),
                    "score": float(hit_data[7]),
                    "bias": float(hit_data[16]),
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
                    "hmmBoundsRaw": hit_data[9],
                    "hmmBounds": hmm_bounds,
                    "evalue": float(hit_data[16]),
                    "score": hit_data[15],
                    "envelopeStart": int(hit_data[13]),
                    "envelopeEnd": int(hit_data[14]),
                    "postProcessed": post_processed,
                    "aliwS": hit_data[6],
                    "..": hit_data[9],
                    "alignment": "",
                    "cigar_alignment": hit_data[17]
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
