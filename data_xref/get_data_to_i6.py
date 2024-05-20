import json
import os

import oracledb


def export_pathways(ipr_uri: str, output_path: str):
    con = oracledb.connect(ipr_uri)
    cur = con.cursor()

    cur.execute(
        """
        SELECT ENTRY_AC, DBCODE, AC, NAME
        FROM INTERPRO.ENTRY2PATHWAY
        """
    )

    pathways = {}
    interpro2pathways = {}
    for entry_acc, dbcode, pathway_id, pathway_name in cur:

        try:
            interpro2pathways[entry_acc].append(pathway_id)
        except KeyError:
            interpro2pathways[entry_acc] = [pathway_id]

        pathways[pathway_id] = [dbcode, pathway_name]

    with open(os.path.join(output_path, "pathways.json"), "wt") as fh:
        json.dump(pathways, fh)

    with open(os.path.join(output_path, "pathways.ipr.json"), "wt") as fh:
        json.dump(interpro2pathways, fh)


def export_goterms(ipr_uri: str, goa_uri: str, output_path: str):
    con = oracledb.connect(ipr_uri)
    cur = con.cursor()

    goa_con = oracledb.connect(goa_uri)
    goa_cur = goa_con.cursor()
    goa_cur.execute(
        """
        SELECT GO_ID, NAME, CATEGORY
        FROM GO.TERMS
        """
    )

    terms = {}
    for go_id, name, category in goa_cur:
        terms[go_id] = [name, category]

    goa_cur.close()
    goa_con.close()

    cur.execute(
        """
        SELECT ENTRY_AC, GO_ID
        FROM INTERPRO.INTERPRO2GO
        WHERE ENTRY_AC IN (
            SELECT ENTRY_AC
            FROM INTERPRO.ENTRY
            WHERE CHECKED = 'Y'
        )
        """
    )

    interpro2go = {}
    for entry_acc, go_id in cur:
        if go_id not in terms:
            continue
        elif entry_acc in interpro2go:
            interpro2go[entry_acc].append(go_id)
        else:
            interpro2go[entry_acc] = [go_id]

    with open(os.path.join(output_path, "goterms.json"), "wt") as fh:
        json.dump(terms, fh)

    with open(os.path.join(output_path, "goterms.ipr.json"), "wt") as fh:
        json.dump(interpro2go, fh)


def export_entries_info(ipr_uri: str, output_path: str):
    con = oracledb.connect(ipr_uri)
    cur = con.cursor()

    cur.execute(
        """
        SELECT M.METHOD_AC, EM.ENTRY_AC, EM.SHORT_NAME, EM.NAME, M.DESCRIPTION, ET.ABBREV
        FROM INTERPRO.METHOD M
        INNER JOIN INTERPRO.CV_DATABASE D
          ON M.DBCODE = D.DBCODE
        INNER JOIN INTERPRO.CV_ENTRY_TYPE ET
          ON M.SIG_TYPE = ET.CODE
        LEFT OUTER JOIN (
            SELECT E.ENTRY_AC, EM.METHOD_AC, E.NAME, E.SHORT_NAME
            FROM INTERPRO.ENTRY E
            INNER JOIN INTERPRO.ENTRY2METHOD EM
              ON E.ENTRY_AC = EM.ENTRY_AC
            WHERE E.CHECKED = 'Y'
        ) EM ON M.METHOD_AC = EM.METHOD_AC
        UNION ALL
        SELECT FM.METHOD_AC, NULL, FM.NAME, NULL, FM.DESCRIPTION, 'Region'
        FROM INTERPRO.FEATURE_METHOD FM
        INNER JOIN INTERPRO.CV_DATABASE D
          ON FM.DBCODE = D.DBCODE
        """
    )

    entries = {}
    for method_ac, entry_ac, short_name, name, description, type in cur:
        entries[method_ac] = [entry_ac, short_name, name, description, type]

    with open(os.path.join(output_path, "entries.ipr.json"), "wt") as fh:
        json.dump(entries, fh)

    cur.close()
    con.close()


if __name__ == '__main__':
    goa_uri = os.environ.get("GOA_URI")
    ippro = os.environ.get("IPR_URI")
    output_path = "./data_xref"

    export_entries_info(ippro, output_path)
    export_pathways(ippro, output_path)
    export_goterms(ippro, goa_uri, output_path)
