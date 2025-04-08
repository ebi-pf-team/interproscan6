import os
import re
import argparse

import oracledb


# Members integrated on InterPro
MEMBER2CODE = {
    "X": "CATH-Gene3D",
    "J": "CDD",
    "Q": "HAMAP",
    "N": "NCBIFAM",
    "V": "PANTHER",
    "H": "Pfam",
    "U": "PIRSF",
    "F": "PRINTS",
    "P": "PROSITE patterns",
    "M": "PROSITE profiles",
    "B": "SFLD",
    "R": "SMART",
    "Y": "SUPERFAMILY"
}


def parse_i6_tsv(i6_output: str):
    interpro2acc = {}
    for member in MEMBER2CODE.values():
        interpro2acc[member] = {}

    with open(i6_output, "rt") as fh:
        lines = fh.readlines()

    for line in lines:
        values = line.split("\t")
        interpro_acc = values[11]
        member_acc = values[4]
        member_db = values[3]
        # Ignore PANTHER subfamilies (won't be integrated)
        if re.fullmatch(r"PTHR\d+:SF\d+", member_acc):
            continue

        try:
            interpro2acc[member_db][interpro_acc].add(member_acc)
        except KeyError:
            interpro2acc[member_db][interpro_acc] = {member_acc}
    return interpro2acc


def swissprot_comparison_i5_i6(ora_url: str, i6_output: str):
    # Get Swiss-Prot <-> Interpro mapping (i6)
    sig2swiss_i6 = parse_i6_tsv(i6_output)

    con = oracledb.connect(ora_url)
    cur = con.cursor()
    info = {}
    cur.execute(
        """
        SELECT M.METHOD_AC, M.DBCODE, E.ENTRY_AC, E.ENTRY_TYPE, E.NAME
        FROM INTERPRO.MATCH M
        INNER JOIN INTERPRO.ENTRY2METHOD EM ON M.METHOD_AC = EM.METHOD_AC
        INNER JOIN INTERPRO.ENTRY E ON EM.ENTRY_AC = E.ENTRY_AC
        INNER JOIN INTERPRO.PROTEIN P ON M.PROTEIN_AC = P.PROTEIN_AC
        WHERE E.CHECKED = 'Y' AND P.DBCODE = 'S'
         and M.DBCODE = 'N'
        """
    )

    sig2swiss_i5 = {}
    for member in MEMBER2CODE.values():
        sig2swiss_i5[member] = {}

    for row in cur.fetchall():
        member_db = MEMBER2CODE[row[1]]
        interpro_acc = row[2]
        member_acc = row[0]
        type = row[3]
        desc = row[4]
        # 'PF06777': ('H', 'IPR010643', 'D', 'Helical and beta-bridge domain')
        info[member_acc] = (member_db, interpro_acc, type, desc)
        try:
            sig2swiss_i5[member_db][row[2]].add(row[0])
        except KeyError:
            sig2swiss_i5[member_db][row[2]] = {row[0]}

    for database in sig2swiss_i5:
        # Compare signatures
        count_lost = 0
        count_gained = 0
        dict_lost = {}
        dict_gained = {}

        for interpro_acc in sig2swiss_i5[database]:
            # Get all members of the signature
            members_i5 = sig2swiss_i5[database][interpro_acc]
            members_i6 = sig2swiss_i6[database].get(interpro_acc, set())

            if members_i5 != members_i6:
                lost = members_i5 - members_i6
                gained = members_i6 - members_i5
                if lost:
                    count_lost += len(lost)
                    try:
                        dict_lost[interpro_acc].update(lost)
                    except KeyError:
                        dict_lost[interpro_acc] = lost
                if gained:
                    count_gained += len(gained)
                    try:
                        dict_gained[interpro_acc].update(lost)
                    except KeyError:
                        dict_gained[interpro_acc] = lost
        print(f"{database}:")
        print(f"  #Lost: {count_lost}")
        if count_lost > 0:
            print(f"    Lost: {dict_lost}")
        print(f"  #Gained: {count_gained}")
        if count_gained > 0:
            print(f"    Gained: {dict_gained}")

    cur.close()
    con.close()


def main():
    parser = argparse.ArgumentParser(prog="IPS_match_regression_test", description="Check presence of matches", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--i6output", type=str, default="tests/data/uniprot_sprot.fasta.tsv", help="TSV with i6 SwissProt result")

    args = parser.parse_args()

    ora_url = os.environ["ORACLE_URL"]
    swissprot_comparison_i5_i6(ora_url, args.i6output)


if __name__ == '__main__':
    main()