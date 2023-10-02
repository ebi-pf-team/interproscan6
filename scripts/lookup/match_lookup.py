import argparse
import json
import xml.etree.ElementTree as ET
from datetime import datetime

import requests

URL_PRECALC_MATCHLOOKUP = "https://www.ebi.ac.uk/interpro/match-lookup/matches"
URL_IS_PRECALC = "https://www.ebi.ac.uk/interpro/match-lookup/isPrecalculated"


def check_precalc(md5: list) -> list:
    all_matches = []
    for i in range(0, len(md5), 4):
        sequences_md5 = ', '.join(md5[i:i + 4])
        checkout = requests.get(f"{URL_IS_PRECALC}?md5={sequences_md5}")
        is_precalc = checkout.text
        if is_precalc:
            match_result = match_lookup(is_precalc)
            all_matches.append(match_result)
    return all_matches


def match_lookup(sequences_md5: str) -> str:
    sequences_input = ('"' + ', '.join(sequences_md5.splitlines()) + '"')
    matches = requests.get(f"{URL_PRECALC_MATCHLOOKUP}?md5={sequences_input}")
    return matches.text


def parse_match(matches: list) -> list[dict]:
    members_info = []
    for match_info in matches:
        tree = ET.fromstring(match_info)

        for match in tree.findall(".//match"):
            protein_md5 = match.find("proteinMD5").text
            for hit in match.findall("hit"):
                hit_data = hit.text.split(",")
                info = {
                    "seq_md5": protein_md5,
                    "analysis": hit_data[0],
                    "signature_acc": hit_data[2],
                    "signature_desc": "",
                    "start": hit_data[4],
                    "stop": hit_data[5],
                    "score": hit_data[16],
                    "date": datetime.today().strftime("%d-%m-%Y"),
                    "match_status": "T",
                    "interpro_annotations_acc": "",
                    "interpro_annotations_desc": "",
                }
                members_info.append(info)

    return members_info


def filter_analysis(members_info: list[dict], applications: list) -> list[dict]:
    if len(applications) > 0:
        filtered_members = []
        applications = list(map(lambda x: x.upper(), applications))
        for member in members_info:
            if member["analysis"].upper() in applications:
                filtered_members.append(member)
        return filtered_members
    else:
        return members_info


def main():
    parser = argparse.ArgumentParser(
        description="Request to precalculated match lookup"
    )
    parser.add_argument(
        "-seq", "--sequences", type=str, help="sequences hash parsed"
    )
    parser.add_argument("-appl", "--applications", nargs="*", help="list of analysis")
    args = parser.parse_args()

    applications = args.applications
    sequences = args.sequences
    not_precalc = []
    json_output = []
    md5 = []
    for seq_id, seq_info in sequences.items():
        md5.append(seq_info[-2])
    md5_upper = [item.upper() for item in md5]
    matches = check_precalc(md5_upper)
    if matches:
        match_parsed = parse_match(matches)
        match_filtered = filter_analysis(match_parsed, applications)
        json_output = json.dumps(match_filtered)
    else:
        not_precalc.append(sequences)
    return json_output, not_precalc


if __name__ == "__main__":
    main()
