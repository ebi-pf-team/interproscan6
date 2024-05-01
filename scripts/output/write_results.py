import json
import os
import sys
from datetime import datetime
import xml.etree.ElementTree as ET


def tsv_output(seq_matches: dict, output_path: str, is_pro: bool):
    def write_to_tsv(
            seq_id, md5, seq_len, member_db, sig_acc,
            ali_from, ali_to, evalue, status, current_date,
            sig_desc="", interpro_acc="-", interpro_desc="-"
    ):
        """
        sig_desc is info from the DB or hmm.out, and added in a later step
        interpro_acc and interpro_desc are added on entries PR
        """
        tsv_file.write((
            f"{seq_id}\t{md5}\t{seq_len}\t{member_db}\t{sig_acc}\t"
            f"{sig_desc}\t{ali_from}\t{ali_to}\t{evalue}\t{status}\t"
            f"{current_date}\t{interpro_acc}\t{interpro_desc}\n"
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
                if match_acc == "signal_peptide":
                    member_db = match["member_db"]
                    sig_acc, status = "Signal Peptide", ""
                    ali_from = match["start"]
                    ali_to = match["end"]
                    evalue = match["pvalue"]
                    write_to_tsv(
                        seq_id, md5, seq_len, member_db,
                        sig_acc, ali_from, ali_to, evalue,
                        status, current_date)
                else:
                    member_db = match["member_db"]
                    sig_acc = match["accession"]
                    status = "T"
                    for location in match["locations"]:
                        evalue = location["evalue"]
                        ali_from = location["start"]
                        ali_to = location["end"]
                        write_to_tsv(
                            seq_id, md5, seq_len, member_db,
                            sig_acc, ali_from, ali_to, evalue,
                            status, current_date)


def json_output(seq_matches: dict, output_path: str, version:str):
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
                        "accession": match_data['accession'],
                        "description": match_data['name'],
                        "signatureLibraryRelease": {
                            "library": match_data['member_db'].upper(),
                            "version": match_data['version']
                        },
                      "entry": match_data['entry']
                    }

                    match = {
                        "signature": signature,
                        "locations": match_data['locations'],
                        "evalue": match_data['evalue'],
                        "score": match_data['score'],
                        "model-ac": match_key
                    }
                matches.append(match)

        result = {
            "sequence": sequence,
            "md5": md5,
            "matches": matches,
            "xref": xrefs
        }
        results.append(result)

    final_data = {"interproscan-version": version, 'results': results}
    with open(json_output, 'w') as json_file:
        json_file.write(json.dumps(final_data, indent=2))

    return final_data


def xml_output(seq_matches: dict, output_path: str, version: str):
    xml_output = os.path.join(output_path + '.xml')
    root = ET.Element("protein-matches", xmlns="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/schemas")
    root.set("interproscan-version", version)

    for seq_id, data in seq_matches.items():
        protein_elem = ET.SubElement(root, "protein")
        sequence_elem = ET.SubElement(protein_elem, "sequence")
        sequence_elem.text = data['sequences'][1]
        sequence_elem.set("md5", data['sequences'][2])
        xref_elem = ET.SubElement(protein_elem, "xref")
        xref_elem.set("id", seq_id)
        xref_elem.set("name", data['sequences'][0])

        matches_elem = ET.SubElement(protein_elem, "matches")
        if 'matches' in data and data['matches']:
            for match_key, match_data in data['matches'].items():
                match_elem = ET.SubElement(matches_elem, "hmmer3-match")
                match_elem.set("evalue", str(match_data['evalue']))
                match_elem.set("score", str(match_data["score"]))

                signature_elem = ET.SubElement(match_elem, "signature")
                signature_elem.set("ac", match_data['accession'])
                signature_elem.set("desc", match_data['description'])
                signature_elem.set("name", match_data['name'])
                if match_data['entry']:
                    entry_elem = ET.SubElement(signature_elem, "entry")
                    entry_elem.set("ac", match_data['entry']['accession'] if match_data['entry']['accession'] else "-")
                    entry_elem.set("desc", match_data['entry']['description'])
                    entry_elem.set("name", match_data['entry']['name'])
                    entry_elem.set("type", match_data['entry']['type'])
                    if match_data['entry']['goXRefs']:
                        for go_xref in match_data['entry']['goXRefs']:
                            go_xref_elem = ET.SubElement(entry_elem, "go-xref")
                            go_xref_elem.set("category", go_xref['category'])
                            go_xref_elem.set("db", go_xref['databaseName'])
                            go_xref_elem.set("id", go_xref['id'])
                            go_xref_elem.set("name", go_xref['name'])
                    if match_data['entry']['pathwayXRefs']:
                        for pathway_xref in match_data['entry']['pathwayXRefs']:
                            pathway_xref_elem = ET.SubElement(entry_elem, "pathway-xref")
                            pathway_xref_elem.set("db", pathway_xref['databaseName'])
                            pathway_xref_elem.set("id", pathway_xref['id'])
                            pathway_xref_elem.set("name", pathway_xref['name'])

                signature_library_elem = ET.SubElement(signature_elem, "signature-library-release")
                signature_library_elem.set("library", match_data['member_db'].upper())
                signature_library_elem.set("version", match_data['version'])
                model_ac_elem = ET.SubElement(match_elem, "model-ac")
                model_ac_elem.text = match_key

                locations_elem = ET.SubElement(match_elem, "locations")
                for location in match_data['locations']:
                    location_elem = ET.SubElement(locations_elem, "hmmer3-location")
                    location_elem.set("env-start", str(location["envelopeStart"]))
                    location_elem.set("env-end", str(location["envelopeEnd"]))
                    location_elem.set("post-processed", str(location["postProcessed"]))
                    location_elem.set("score", str(location["score"]))
                    location_elem.set("evalue", str(location["evalue"]))
                    location_elem.set("hmm-start", str(location["hmmStart"]))
                    location_elem.set("hmm-end", str(location["hmmEnd"]))
                    location_elem.set("hmm-length", str(location["hmmLength"]))
                    location_elem.set("hmm-bounds", str(location["hmmBounds"]))
                    location_elem.set("start", str(location["start"]))
                    location_elem.set("end", str(location["end"]))
                    location_elem.set("representative", str(location["representative"]))
                    location_frags_elem = ET.SubElement(location_elem, "location-fragments")
                    if 'sites' in location:
                        location_frag_elem = ET.SubElement(location_frags_elem, "hmmer3-location-fragment")
                        location_frag_elem.set("start", str(location["sites"]["start"]))
                        location_frag_elem.set("end", str(location["sites"]["end"]))
                        location_frag_elem.set("dc-status", "")

    tree = ET.ElementTree(root)
    tree.write(xml_output, encoding="utf-8", xml_declaration=True)


def write_results(sequences_path: str, matches_path: str, output_format: str, output_path: str, version: str):
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

    if "TSV" in output_format:
        tsv_output(seq_matches, output_path, False)
    if "TSV-PRO" in output_format:
        tsv_output(seq_matches, output_path, True)
    if "JSON" in output_format:
        json_output(seq_matches, output_path, version)
    if "XML" in output_format:
        xml_output(seq_matches, output_path, version)


def main():
    args = sys.argv[1:]

    sequences = args[0]
    matches = args[1]
    formats = args[2]
    output_path = args[3]
    version = args[4]

    write_results(sequences, matches, formats, output_path, version)


if __name__ == "__main__":
    main()
