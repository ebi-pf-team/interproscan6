import argparse
import ast
import json
import xml.etree.ElementTree as ET
from datetime import datetime

import requests


def match_lookup(sequences_md5: str, url: str) -> str:
    with open(sequences_md5, 'r') as md5_data:
        md5_info = md5_data.read()
    checked_seq_md5 = ast.literal_eval(md5_info)
    matches = checked_seq_md5["matches"]
    url_input = ','.join(matches)

    matches = requests.get(f"{url}?md5={url_input}")
    return matches.text


def parse_match(matches: str) -> list[dict]:
    members_info = []
    # for match_info in matches:
    tree = ET.fromstring(matches)

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
        "-checked", "--checked_lookup", type=str, help="dict with md5 lookup matches checked"
    )
    parser.add_argument("-appl", "--applications", type=str, help="list of analysis")
    parser.add_argument("-url", "--url", type=str, help="url to get sequences match lookup")
    args = parser.parse_args()

    applications = args.applications[1: -1].split(', ')
    match_results = match_lookup(args.checked_lookup, args.url)
    match_parsed = parse_match(match_results)
    match_filtered = filter_analysis(match_parsed, applications)
    json_output = json.dumps(match_filtered)
    print(json_output)


if __name__ == "__main__":
    main()
