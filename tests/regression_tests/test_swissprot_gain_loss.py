import os
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
    entry2proteins_i5 = {}

    con = oracledb.connect(ora_url)
    cur = con.cursor()

    cur.execute(
        """
        SELECT DISTINCT E.ENTRY_AC, M.PROTEIN_AC
        FROM INTERPRO.MATCH M
        INNER JOIN INTERPRO.ENTRY2METHOD EM ON M.METHOD_AC = EM.METHOD_AC
        INNER JOIN INTERPRO.ENTRY E ON EM.ENTRY_AC = E.ENTRY_AC
        INNER JOIN INTERPRO.PROTEIN P ON M.PROTEIN_AC = P.PROTEIN_AC
        WHERE E.CHECKED = 'Y' AND P.DBCODE = 'S'
        """
    )

    for entry, protein in cur.fetchall():
        try:
            entry2proteins_i5[entry].append(protein)
        except KeyError:
            entry2proteins_i5[entry] = [protein]

    cur.close()
    con.close()

    return entry2proteins_i5


def get_i6_swissprot_mapping(i6_result: dict):
    entry2proteins = {}
    for result in i6_result["results"]:
        for match in result["matches"]:
            entry = match["signature"]["entry"]
            if not entry:
                continue
            else:
                entry_acc = match["signature"]["entry"]["accession"]
                for protein in result["xref"]:
                    protein_id = protein["id"].split("|")[1]
                    try:
                        entry2proteins[entry_acc].append(protein_id)
                    except KeyError:
                        entry2proteins[entry_acc] = [protein_id]
    return entry2proteins


def swissprot_comparison(entry2proteins_i5: dict, entry2proteins_i6: dict):
    count_total_lost = 0
    count_total_gained = 0
    for entry in entry2proteins_i5:
        members_i6 = set(entry2proteins_i6.get(entry, set()))
        members_i5 = set(entry2proteins_i5.get(entry, set()))
        if members_i5 != members_i6:
            lost = members_i5 - members_i6
            gained = members_i6 - members_i5
            if lost or gained:
                count_total_lost += len(lost)
                count_total_gained += len(gained)
                print(f"{entry} | #Total: {len(members_i5)} | #Lost: {len(lost)} | #Gained: {len(gained)} | {lost} | {gained}")
    print(f"Total lost: {count_total_lost} | Total gained: {count_total_gained}")


def main():
    parser = argparse.ArgumentParser(prog="IPS_swissprot_test", description="Check gain/lost of swissprot between i5 and i6", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--i6output", type=str, help="i6 SwissProt JSON result")

    args = parser.parse_args()
    ora_url = os.environ["ORACLE_URL"]

    entry2proteins_i5 = get_i5_swissprot_mapping(ora_url)

    with open(args.i6output, "r") as f:
        entry2proteins_i6 = get_i6_swissprot_mapping(json.load(f))

    swissprot_comparison(entry2proteins_i5, entry2proteins_i6)


if __name__ == '__main__':
    main()