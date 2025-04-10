import os
import re
import argparse
import json

import oracledb


# Members integrated on InterPro
PREFIX2MEMBERDB = {
    "G3DSA": "CATH-Gene3D",
    "cd": "CDD",
    "MF": "HAMAP",
    "NF": "NCBIFAM",
    "TIGR": "NCBIFAM (TIGRFAM)",
    "PTHR": "PANTHER",
    "PF": "Pfam",
    "PIRSF": "PIRSF",
    "PR": "PRINTS",
    "PS": "PROSITE",
    "SFLD": "SFLD",
    "SM": "SMART",
    "SSF": "SUPERFAMILY"
}


def get_i5_swissprot_mapping(ora_url: str):
    entry2sigs_i5 = {}

    con = oracledb.connect(ora_url)
    cur = con.cursor()

    cur.execute(
        """
        SELECT E.ENTRY_AC, M.METHOD_AC
        FROM INTERPRO.MATCH M
        INNER JOIN INTERPRO.ENTRY2METHOD EM ON M.METHOD_AC = EM.METHOD_AC
        INNER JOIN INTERPRO.ENTRY E ON EM.ENTRY_AC = E.ENTRY_AC
        INNER JOIN INTERPRO.PROTEIN P ON M.PROTEIN_AC = P.PROTEIN_AC
        WHERE E.CHECKED = 'Y' AND P.DBCODE = 'S'
        """
    )

    for entry, signature in cur.fetchall():
        try:
            entry2sigs_i5[entry].add(signature)
        except KeyError:
            entry2sigs_i5[entry] = {signature}

    cur.close()
    con.close()

    return entry2sigs_i5


def get_i6_swissprot_mapping(i6_result: dict):
    entry2sigs = {}

    for protein_dict in i6_result["results"]:
        for match in protein_dict["matches"]:
            try:
                entry = match["signature"]["entry"]
                if not entry:
                    continue
                else:
                    entry_acc = match["signature"]["entry"]["accession"]
            except TypeError:
                print(f"KeyError: {match['signature']}")
            signature = match["signature"]["accession"]

            # Ignore PANTHER subfamilies (won't be integrated)
            if re.fullmatch(r"PTHR\d+:SF\d+", signature):
                continue

            try:
                entry2sigs[entry_acc].add(signature)
            except KeyError:
                entry2sigs[entry_acc] = {signature}

    return entry2sigs


def swissprot_comparison(entry2sigs_i5: dict, entry2sigs_i6: dict):
    count_total_lost = 0
    count_total_gained = 0
    for entry in entry2sigs_i5:
        members_i6 = entry2sigs_i6.get(entry, set())
        members_i5 = entry2sigs_i5.get(entry, set())
        if members_i5 != members_i6:
            lost = members_i5 - members_i6
            gained = members_i6 - members_i5
            if lost or gained:
                count_total_lost += len(lost)
                count_total_gained += len(gained)
                print(f"Entry: {entry} | #Lost: {len(lost)} | #Gained: {len(gained)} | Lost: {lost} | Gained: {gained}")
    print(f"Total lost: {count_total_lost} | Total gained: {count_total_gained}")


def main():
    parser = argparse.ArgumentParser(prog="IPS_swissprot_test", description="Check gain/lost of swissprot between i5 and i6", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--i6output", type=str, help="i6 SwissProt JSON result")
    parser.add_argument("--bymembers", help="shows gain/lost by memberDB", action="store_true")

    args = parser.parse_args()
    ora_url = os.environ["ORACLE_URL"]

    with open(args.i6output, "r") as fh:
        entry2sigs_i6 = get_i6_swissprot_mapping(json.load(fh))

    entry2sigs_i5 = get_i5_swissprot_mapping(ora_url)

    if args.bymembers:
        # Compare by memberDB
        for prefix in PREFIX2MEMBERDB:
            entry2sig_i5_by_member = {}
            entry2sig_i6_by_member = {}
            for entry, signatures in entry2sigs_i5.items():
                entry2sig_i5_by_member[entry] = {signature for signature in signatures if signature.startswith(prefix)}
            for entry, signatures in entry2sigs_i6.items():
                entry2sig_i6_by_member[entry] = {signature for signature in signatures if signature.startswith(prefix)}
            print(f"{PREFIX2MEMBERDB[prefix]}:")
            swissprot_comparison(entry2sig_i5_by_member, entry2sig_i6_by_member)
    else:
        swissprot_comparison(entry2sigs_i5, entry2sigs_i6)


if __name__ == '__main__':
    main()