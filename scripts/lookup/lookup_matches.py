import argparse
import ast
import json
import xml.etree.ElementTree as ET

import requests


def match_lookup(matches_checked: list, url: str) -> str:
    url_input = ','.join(matches_checked)
    matches = requests.get(f"{url}?md5={url_input}")
    return matches.text


def parse_match(matches: str, applications: list, md52seq_id: dict, match_parsed: dict) -> dict:
    member_matches = []
    tree = ET.fromstring(matches)

    for match in tree.findall(".//match"):
        match_id = match.find("matchId").text
        protein_md5 = match.find("proteinMD5").text
        hits = []

        for hit in match.findall("hit"):
            hit_data = hit.text.split(',')
            hit_appl = hit_data[0]
            if hit_appl in applications:
                domain = {
                    "application": hit_appl,
                    "version": hit_data[1],
                    "accession": hit_data[2],
                    "accession2": hit_data[3],
                    "ali_from": hit_data[4],
                    "ali_to": hit_data[5],
                    "aliwS": hit_data[6],
                    "score": hit_data[7],
                    "bias": hit_data[8],
                    "..": hit_data[9],
                    "hmm_from": hit_data[10],
                    "hmm_to": hit_data[11],
                    "qlen": hit_data[12],
                    "env_from": hit_data[13],
                    "env_to": hit_data[14],
                    "score_hit": hit_data[15],
                    "e_value": hit_data[16],
                    "cigar_alignment": hit_data[17],
                }
                hits.append(domain)

        if hits:
            match_dict = {
                "match_id": match_id,
                "md5": protein_md5,
                "domains": hits
            }

            member_matches.append(match_dict)
            seq_id = md52seq_id[protein_md5]
            try:
                match_parsed[seq_id].append(member_matches)
            except:
                match_parsed[seq_id] = member_matches

    return match_parsed


def main():
    parser = argparse.ArgumentParser(
        description="Parse to match lookup"
    )
    parser.add_argument(
        "-checked", "--checked_lookup", type=str, help="dict with md5 lookup matches checked"
    )
    parser.add_argument("-appl", "--applications", type=str, help="list of analysis")
    parser.add_argument("-url", "--url", type=str, help="url to get sequences match lookup")
    args = parser.parse_args()

    applications = list(map(lambda x: x.upper(), args.applications[1: -1].split(', ')))
    match_parsed = {}

    with open(args.checked_lookup, 'r') as md5_data:
        checked_data = json.load(md5_data)
    matches = checked_data["matches"]
    seq_info = checked_data["sequences_info"]

    md52seq_id = {}
    for seq_id, match in seq_info.items():
        md52seq_id[match[-2]] = seq_id

    match_results = match_lookup(matches, args.url)
    match_parsed = parse_match(match_results, applications, md52seq_id, match_parsed)

    print(json.dumps(match_parsed))


if __name__ == "__main__":
    main()
