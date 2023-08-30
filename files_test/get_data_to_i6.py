import json
import os

import oracledb


def _export_pathways(cur: oracledb.Cursor, output_path: str):
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


def _export_go_terms(cur: oracledb.Cursor, goa_uri: str, output_path: str):
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


def export_for_interproscan(ipr_uri: str, goa_uri: str, outdir: str):
    con = oracledb.connect(ipr_uri)
    cur = con.cursor()

    _export_pathways(cur, outdir)
    _export_go_terms(cur, goa_uri, outdir)

    cur.close()
    con.close()


def export_interpro_info(ipr_uri: str, output_path: str):
    con = oracledb.connect(ipr_uri)
    cur = con.cursor()

    cur.execute(
        """SELECT DISTINCT M.METHOD_AC, NAME, DESCRIPTION, ENTRY_AC
            FROM INTERPRO.METHOD M
            LEFT OUTER JOIN (
                SELECT E.ENTRY_AC, EM.METHOD_AC
                FROM INTERPRO.ENTRY E
                INNER JOIN INTERPRO.ENTRY2METHOD EM
                  ON E.ENTRY_AC = EM.ENTRY_AC
                WHERE E.CHECKED = 'Y'
            ) EM ON M.METHOD_AC = EM.METHOD_AC
        """
    )

    entries = {}
    for method_ac, name, description, entry_ac in cur:
        entries[method_ac] = [name, description, entry_ac]

    with open(os.path.join(output_path, "ipr_entries.json"), "wt") as fh:
        json.dump(entries, fh)

    cur.close()
    con.close()


if __name__ == '__main__':
    goa_uri = ""
    ippro = ""
    output_path = "./i6data"
    export_for_interproscan(ippro, goa_uri, output_path)
    export_interpro_info(ippro, output_path)
