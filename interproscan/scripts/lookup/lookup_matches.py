import json
import logging
import sys
import urllib.request
import xml.etree.ElementTree as ET

from retry_conn_decorator import lookup_retry_decorator


@lookup_retry_decorator
def match_lookup(matches_checked: list, url: str, **kwargs) -> str:
    url_input = ','.join(matches_checked)
    matches = urllib.request.urlopen(f"{url}?md5={url_input}")
    print(f"{url}?md5={url_input}")
    return matches.read().decode('utf-8')


def parse_match(match_data: str, applications: list, md52seq_id: dict) -> dict:
    tree = ET.fromstring(match_data)
    matches = {}
    hmm_bound_pattern = {"[]": "COMPLETE", "[.": "N_TERMINAL_COMPLETE", ".]": "C_TERMINAL_COMPLETE", "..": "INCOMPLETE"}
    """
    Structure of hit data:
    1. Member db version
    2. Hit Accession
    3. Model accession
    4. Location start
    5. Location end
    6. Location fragment
    ...17
    """

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

                if signature["member_db"].upper() == "PANTHER":
                    signature["node_id"] = hit_data[17]

                if signature["member_db"].upper() == "MOBIDBLITE":
                    signature["sequence-feature"] = hit_data[17]

                location = {
                    "start": int(hit_data[4]),
                    "end": int(hit_data[5]),
                    "representative": "",
                    "hmmStart": int(hit_data[10]),
                    "hmmEnd": int(hit_data[11]),
                    "hmmLength": int(hit_data[12]),
                    "hmmBoundsRaw": hit_data[9],
                    "hmmBounds": hmm_bounds,
                    "score": hit_data[15],
                    "envelopeStart": int(hit_data[13]),
                    "envelopeEnd": int(hit_data[14]),
                    "postProcessed": post_processed,
                    "locationFragment": hit_data[6],
                    # misc either , 0 or HmmBounds raw
                    "misc": hit_data[9],
                    "alignment": "",
                    "cigar_alignment": hit_data[17]
                }

                if hit_appl != "PRINTS":
                    location["evalue"] = float(hit_data[16])
                else:
                    location["pvalue"] = float(hit_data[16])
                    # prints mls stores motif number at same index as hmm length
                    # set hmm length to 0
                    location["hmmLength"] = int(0)
                    location["motifNumber"] = hit_data[12]
                    signature["graphscan"] = hit_data[17]

                if target_key not in matches:
                    matches[target_key] = {}

                if accession not in matches[target_key]:
                    matches[target_key][accession] = signature
                    matches[target_key][accession]["locations"] = [location]
                else:
                    matches[target_key][accession]["locations"].append(location)

    return matches


def main():
    """CL input:
    0. Str repr of path to a JSON file containing the MD5 hashed sequences
    1. Str repr of list of selected application
    2. Url for MLS
    3. Num of times to retry the connection
    4. Str repr of path for the output file
    """
    args = sys.argv[1:]
    checked_lookup = args[0]
    applications = args[1].split(',')
    url = args[2]
    retries = int(args[3])

    applications = list(map(lambda x: x.upper().replace('MOBIDB', 'MOBIDB_LITE'), applications))

    with open(checked_lookup, 'r') as md5_data:
        checked_data = json.load(md5_data)
    matches = checked_data["matches"]
    seq_info = checked_data["sequences_info"]

    md52seq_id = {}
    for seq_id, match in seq_info.items():
        md52seq_id[match['md5']] = seq_id

    match_results, err = match_lookup(matches, url, retries=retries)

    if err:
        logging.error(err)

    if match_results:
        match_parsed = parse_match(match_results, applications, md52seq_id)
        if len(match_parsed) > 0:
            with open(args[4], "w") as fh:
                json.dump(match_parsed, fh)


if __name__ == "__main__":
    main()
