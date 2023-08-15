import argparse
import hashlib
import json
import xml.etree.ElementTree as ET
from datetime import datetime

import requests

URL_PRECALC_MATCHLOOKUP = "https://www.ebi.ac.uk/interpro/match-lookup/matches"


def get_sequences(fasta_file: str) -> dict:
    sequences = {}
    with open(fasta_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                seq_id = line.split(" ")[0][1:]
                sequences[seq_id] = ""
            else:
                sequences[seq_id] += line
    return sequences


def match_lookup(sequence: str) -> str:
    sequence_md5 = hashlib.md5(sequence.encode()).hexdigest().upper()
    matches = requests.get(f"{URL_PRECALC_MATCHLOOKUP}?md5={sequence_md5}")
    return matches.text


def parse_match(id: str, sequence: str, matches: str) -> list[dict]:
    members_info = []
    tree = ET.fromstring(matches)

    for match in tree.findall(".//match"):
        protein_md5 = match.find("proteinMD5").text
        for hit in match.findall("hit"):
            hit_data = hit.text.split(",")
            info = {
                "protein_acc": id,
                "seq_md5": protein_md5,
                "seq_len": len(sequence),
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
    parser.add_argument("-fasta", "--fastafile", type=str, help="list of analysis")
    parser.add_argument("-appl", "--applications", nargs="*", help="list of analysis")
    args = parser.parse_args()

    applications = args.applications

    sequences = get_sequences(args.fastafile)
    for id, seq in sequences.items():
        matches = match_lookup(seq)
        match_parsed = parse_match(id, seq, matches)
        match_filtered = filter_analysis(match_parsed, applications)
        print(json.dumps(match_filtered))


if __name__ == "__main__":
    main()
