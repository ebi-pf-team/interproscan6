import json
import os
import sys
from datetime import datetime


def tsv_output(seq_matches: dict, output_path: str, is_pro: bool):
    def write_to_tsv(
            seq_id, md5, seq_len, member_db, sig_acc,
            sig_desc, ali_from, ali_to, evalue, status,
            current_date, interpro_acc, interpro_name, xrefs
    ):
        tsv_file.write((
            f"{seq_id}\t{md5}\t{seq_len}\t{member_db}\t{sig_acc}\t"
            f"{sig_desc}\t{ali_from}\t{ali_to}\t{evalue}\t{status}\t"
            f"{current_date}\t{interpro_acc}\t{interpro_name}\t"
            f"{xrefs}\n"
        ))

    tsv_output = os.path.join(output_path + '.tsv')
    if is_pro:
        tsv_output = tsv_output + "-pro"

    with open(tsv_output, 'w') as tsv_file:
        current_date = datetime.now().strftime('%d-%m-%Y')
        for seq_target, info in seq_matches.items():
            sequence_data = info['sequences']
            matches = info["matches"]
            seq_id = seq_target
            md5 = sequence_data[2]
            seq_len = sequence_data[3]

            for match_acc, match in matches.items():
                match_db = match["member_db"]
                entry_acc = match["entry"]["accession"]
                entry_name = match["entry"]["name"]
                entry_desc = match["entry"]["description"]
                goterms = []
                pathways = []
                for go_info in match["entry"]["goXRefs"]:
                    goterms.append(go_info["id"])
                for pwy_info in match["entry"]["pathwayXRefs"]:
                    pathways.append(pwy_info["id"])
                xrefs = f"{'|'.join(goterms)}\t{'|'.join(pathways)}"

                if match_acc == "signal_peptide":
                    sig_acc, status = "Signal Peptide", ""
                    ali_from = match["start"]
                    ali_to = match["end"]
                    evalue = match["pvalue"]
                else:
                    sig_acc = match["accession"]
                    status = "T"
                    for location in match["locations"]:
                        evalue = location["evalue"]
                        ali_from = location["start"]
                        ali_to = location["end"]
                write_to_tsv(
                    seq_id, md5, seq_len, match_db,
                    sig_acc, entry_desc, ali_from, ali_to,
                    evalue, status, current_date, entry_acc,
                    entry_name, xrefs)


def json_output(seq_matches: dict, output_path: str):
    json_output = os.path.join(output_path + '.json')
    results = []

    for seq_id, data in seq_matches.items():
        sequence = data['sequences'][1]
        md5 = data['sequences'][2]
        matches = []
        xrefs = {
            "name": data['sequences'][0],
            "id": seq_id
        }
        if 'matches' in data and data['matches']:
            for match_key, match_data in data['matches'].items():  # match_key == sig_Acc
                if match_key == "signal_peptide":
                    match = {
                        "signature": match_key,
                        "SignalP_release": match_data["signalp_version"],
                        "start": match_data["start"],
                        "end": match_data["end"],
                        "pvalue": match_data["pvalue"],
                    }
                else:
                    signature = {
                        "accession": match_data['accession'].split(":")[0],  # drop subfamily
                        "name": match_data['name'],
                        "description": match_data["entry"]["description"],
                        "signatureLibraryRelease": {
                            "library": match_data['member_db'].upper(),
                            "version": match_data['version']
                        },
                        "entry": match_data['entry'] if match_data['entry']['accession'] != "-" else None
                    }

                    match = {
                        "signature": signature,
                        "locations": match_data['locations'],
                    }

                    if match_data['member_db'].upper() == "CDD":
                        match["evalue"] = match_data['locations'][0]["evalue"]
                        match["score"] = match_data['locations'][0]["score"]

                    else:
                        match["evalue"] = match_data['evalue']
                        match["score"] = match_data['score']

                    match["model-ac"] = match_data['model-ac']

                    if match_data['member_db'].upper() == "PANTHER":
                        # get protein class and graftpoint for Panther
                        match['proteinClass'] = match_data['proteinClass']
                        match['graftPoint'] = match_data['graftPoint']

                matches.append(match)

        result = {
            "sequence": sequence,
            "md5": md5,
            "matches": matches,
            "xref": xrefs
        }
        results.append(result)

    final_data = {"interproscan-version": "6.0.0", 'results': results}
    with open(json_output, 'w') as json_file:
        json_file.write(json.dumps(final_data, indent=2))


def write_results(sequences_path: str, matches_path: str, output_format: str, output_path: str):
    output_format = output_format.upper()
    seq_matches = {}

    all_sequences = {}
    with open(matches_path, 'r') as match_data:
        all_matches = json.load(match_data)
    with open(sequences_path, 'r') as seq_data:
        for line in seq_data:
            sequence = json.loads(line)
            all_sequences.update(sequence)

    for key in all_sequences:
        if key in all_matches:
            seq_matches[key] = {
                'sequences': all_sequences[key],
                'matches': all_matches[key]
            }
        else:
            seq_matches[key] = {
                'sequences': all_sequences[key],
                'matches': {}
            }

    print(json.dumps(seq_matches, indent=4))

    if "TSV" in output_format:
        tsv_output(seq_matches, output_path, False)
    if "TSV-PRO" in output_format:
        tsv_output(seq_matches, output_path, True)
    if "JSON" in output_format:
        json_output(seq_matches, output_path)


def main():
    args = sys.argv[1:]

    sequences = args[0]
    matches = args[1]
    formats = args[2]
    output_path = args[3]

    write_results(sequences, matches, formats, output_path)


if __name__ == "__main__":
    main()
