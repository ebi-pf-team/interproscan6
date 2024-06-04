import json
import os
import sys
from datetime import datetime
import xml.etree.ElementTree as ET


MATCH_ELEMENT = {
    'SIGNALP': 'signal-peptide',
    'CDD': 'cdd-domain',
    'ANTIFAM': 'hmmer3-match',
    'GENE3D': 'hmmer3-match',
    'FUNFAM': 'hmmer3-match',
    'NCBIFAM': 'hmmer3-match',
    'PANTHER': 'hmmer3-match',
    'SFLD': 'hmmer3-match',
}


def tsv_output(seq_matches: dict, output_path: str):
    def write_to_tsv(
            seq_id, md5, seq_len, member_db, sig_acc,
            sig_desc, ali_from, ali_to, evalue, status,
            current_date, interpro_acc, interpro_name, xrefs):
        tsv_file.write((
            f"{seq_id}\t{md5}\t{seq_len}\t{member_db}\t{sig_acc}\t"
            f"{sig_desc}\t{ali_from}\t{ali_to}\t{evalue}\t{status}\t"
            f"{current_date}\t{interpro_acc}\t{interpro_name}\t{xrefs}\n"
        ))

    tsv_output = os.path.join(output_path + '.tsv')

    with open(tsv_output, 'w') as tsv_file:
        current_date = datetime.now().strftime('%d-%m-%Y')
        for seq_target, info in seq_matches.items():
            sequence_data = info['sequences']
            matches = info["matches"]
            seq_id = seq_target
            md5 = sequence_data[2]
            seq_len = sequence_data[3]

            for match_acc, match in matches.items():
                entry_acc, entry_name, entry_desc = "-", "-", "-"
                goterms, pathways = [], []
                if match["entry"]:
                    entry_acc = match["entry"]["accession"]
                    entry_name = match["entry"]["short_name"]
                    entry_desc = match["entry"]["name"]
                    for go_info in match["entry"]["goXRefs"]:
                        goterms.append(go_info["id"])
                    for pwy_info in match["entry"]["pathwayXRefs"]:
                        pathways.append(pwy_info["id"])
                match_db = match["member_db"]
                xrefs = f"{'|'.join(goterms)}\t{'|'.join(pathways)}"

                for location in match["locations"]:
                    if match_acc == "signal_peptide":
                        sig_acc, status = "Signal Peptide", ""
                        ali_from = match["locations"][0]["start"]
                        ali_to = match["locations"][0]["end"]
                        evalue = match["locations"][0]["pvalue"]
                    else:
                        sig_acc = match["accession"]
                        status = "T"
                        evalue = location["evalue"]
                        ali_from = location["start"]
                        ali_to = location["end"]

                    write_to_tsv(
                        seq_id, md5, seq_len, match_db,
                        sig_acc, entry_desc, ali_from, ali_to,
                        evalue, status, current_date, entry_acc,
                        entry_name, xrefs)


def tsv_pro_output(seq_matches: dict, output_path: str):
    def write_to_tsv(
            member_db, version_major, version_minor, seq_id,
            sig_acc, model_ac, ali_from, ali_to, fragment,
            score, evalue, raw_hmm_bound, hmm_start, hmm_end,
            hmm_length, env_start, env_end, location_score,
            location_evalue, cigar_alignment):
        tsv_file.write((
            f"{member_db}\t{version_major}\t{version_minor}\t{seq_id}\t"
            f"{sig_acc}\t{model_ac}\t{ali_from}\t{ali_to}\t{fragment}\t"
            f"{score}\t{evalue}\t{raw_hmm_bound}\t{hmm_start}\t{hmm_end}\t"
            f"{hmm_length}\t{env_start}\t{env_end}\t{location_score}\t"
            f"{location_evalue}\t{cigar_alignment}\n"
        ))

    tsv_output = os.path.join(output_path + '.tsv-pro')

    with open(tsv_output, 'w') as tsv_file:
        for seq_target, info in seq_matches.items():
            matches = info["matches"]
            seq_id = seq_target

            for match_acc, match in matches.items():
                member_db = match["member_db"]
                try:
                    version_major, version_minor = match['version'].split('.')
                except ValueError:
                    version_major = match['version']
                    version_minor = "0"
                try:  # cdd does not have evalue and score on this level
                    evalue = match["evalue"]
                    score = match["score"]
                except KeyError:
                    evalue = "-"
                    score = "-"

                if 'model-ac' in match:
                    model_ac = match['model-ac']
                elif member_db.upper() in ["SIGNALP"]:  # will probably apply to TMHMM and Phobius when added
                    model_ac = "-"
                else:
                    model_ac = match['accession']

                for location in match["locations"]:
                    if member_db.upper() == "CDD":
                        hmm_start = location["start"]
                        hmm_end = location["end"]
                        hmm_length = int(hmm_end) - int(hmm_start)
                        location_score = location["score"]
                        env_end, env_start = "-", "-"
                    elif member_db.upper() == "SIGNALP":
                        hmm_start = location["start"]
                        hmm_end = location["end"]
                        hmm_length = int(hmm_end) - int(hmm_start)
                        location_score = location["pvalue"]
                        env_end, env_start = "-", "-"
                    else:
                        hmm_start = location["hmmStart"]
                        hmm_end = location["hmmEnd"]
                        hmm_length = location["hmmLength"]
                        location_score = location["score"]
                        env_end = location["envelopeEnd"]
                        env_start = location["envelopeStart"]
                    try:
                        fragment = location["fragment"]
                    except KeyError:
                        fragment = f"{location['start']}-{location['end']}-S"
                    try:
                        raw_hmm_bound = location["rawHmmBounds"]
                    except KeyError:
                        raw_hmm_bound = ""  # lookup match does not have rawHmmBounds

                    if match_acc == "signal_peptide":
                        sig_acc, status = "Signal Peptide", ""
                        ali_from = match["locations"][0]["start"]
                        ali_to = match["locations"][0]["end"]
                        location_evalue = match["locations"][0]["pvalue"]
                    else:
                        sig_acc = match["accession"]
                        evalue = location["evalue"]
                        ali_from = location["start"]
                        ali_to = location["end"]
                        location_evalue = location["evalue"]
                    cigar_alignment = ""
                    try:
                        cigar_alignment = location["cigar_alignment"]
                    except KeyError:
                        pass  # some members may not have cigar alignment

                    write_to_tsv(
                        member_db, version_major, version_minor, seq_id,
                        sig_acc, model_ac, ali_from, ali_to, fragment,
                        score, evalue, raw_hmm_bound, hmm_start, hmm_end,
                        hmm_length, env_start, env_end, location_score,
                        location_evalue, cigar_alignment)


def json_output(seq_matches: dict, output_path: str, version: str):
    json_output = os.path.join(output_path + '.json')
    results = []
    boolean_map = {"true": True, "false": False}
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
                        "SignalP_release": match_data["version"],
                        "start": match_data["locations"][0]["start"],
                        "end": match_data["locations"][0]["end"],
                        "pvalue": match_data["locations"][0]["pvalue"],
                    }
                else:
                    try:
                        description = match_data['description']
                    except KeyError:
                        description = "-"
                    entry = None
                    if match_data['entry']['accession'] != "-":
                        if match_data['member_db'].upper() == "PANTHER":
                            entry = {
                                "name": match_data['entry']['description']
                            }
                        else:
                            description = match_data['entry']['description']
                            entry = {
                                "accession": match_data['entry']['accession'],
                                "name": match_data['entry']['short_name'],
                                "description": match_data['entry']['name'],
                                "type": match_data['entry']['type'].upper()
                            }
                        try:
                            entry["goXRefs"] = match_data['entry']['goXRefs']
                        except KeyError:
                            if match_data['member_db'].upper() == "PANTHER":
                                entry["goXRefs"] = []
                            else:
                                pass
                        try:
                            entry["pathwayXRefs"] = match_data['entry']['pathwayXRefs']
                        except KeyError:
                            pass
                    if match_data['member_db'].upper() == "CDD":
                        description = match_data['name']
                    signature = {
                        "accession": match_data['accession'].split(":")[0],  # drop subfamily
                        "name": match_data['name'],
                        "description": description,
                        "signatureLibraryRelease": {
                            "library": match_data['member_db'].upper(),
                            "version": match_data['version']
                        },
                        "entry": entry
                    }
                    locations = []
                    for location in match_data['locations']:
                        info = {
                            "start": int(location["start"]),
                            "end": int(location["end"]),
                            "representative": boolean_map.get(location["representative"].lower(), False),
                            "evalue": float(location["evalue"]),
                            "score": float(location["score"]),
                        }
                        if match_data['member_db'].upper() == "CDD":
                            pass
                        else:
                            info["hmmStart"] = int(location["hmmStart"])
                            info["hmmEnd"] = int(location["hmmEnd"])
                            info["hmmLength"] = int(location["hmmLength"])
                            info["hmmBounds"] = location["hmmBounds"]
                            info["envelopeStart"] = int(location["envelopeStart"])
                            info["envelopeEnd"] = int(location["envelopeEnd"])
                            info["postProcessed"] = boolean_map.get(location["postProcessed"].lower())
                        if match_data['member_db'].upper() in ["SFLD", "CDD"]:
                            info["sites"] = location["sites"]
                        try:
                            info["location-fragments"] = location["location-fragments"]
                        except KeyError:
                            single_location = {
                                "start": int(location["start"]),
                                "end": int(location["end"]),
                                "dc-status": "CONTINUOUS"
                            }
                            try:
                                info["location-fragments"].append(single_location)
                            except KeyError:
                                info["location-fragments"] = [single_location]
                        locations.append(info)
                    match = {
                        "signature": signature,
                        "locations": locations
                    }

                    if match_data['member_db'].upper() != "CDD":
                        match["evalue"] = float(match_data['evalue'])
                        match["score"] = float(match_data['score'])

                    if 'model-ac' in match_data:
                        match["model-ac"] = match_data['model-ac']
                    else:
                        match["model-ac"] = match_data['accession']

                    if match_data['member_db'].upper() == "PANTHER":
                        # get protein class and graftpoint for Panther
                        match['accession'] = match_data['accession']
                        try:
                            match['proteinClass'] = match_data['proteinClass']
                            match['graftPoint'] = match_data['graftPoint']
                        except KeyError:
                            pass

                matches.append(match)

        result = {
            "sequence": sequence,
            "md5": md5,
            "matches": matches,
            "xref": [xrefs]
        }
        results.append(result)

    final_data = {"interproscan-version": version, 'results': results}
    with open(json_output, 'w') as json_file:
        json_file.write(json.dumps(final_data, indent=2))

    return final_data


def xml_output(seq_matches: dict, output_path: str, version: str):
    xml_output = os.path.join(output_path + '.xml')
    root = ET.Element("protein-matches", xmlns="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6/schemas")
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
                match_elem = ET.SubElement(matches_elem, MATCH_ELEMENT[match_data['member_db'].upper()])

                try:
                    match_elem.set("evalue", str(match_data['evalue']).upper())
                    match_elem.set("score", str(match_data["score"]))
                except KeyError:
                    pass  # some members may not have evalue or score on this level (e.g. cdd)

                signature_elem = ET.SubElement(match_elem, "signature")
                if match_data['member_db'].upper() not in ['SIGNALP']:  # member db that don't have sigs, so no accs etc.
                    signature_elem.set("ac", match_data['accession'])
                    signature_elem.set("desc", match_data['name'])
                    signature_elem.set("name", match_data['name'])

                if match_data['entry']:
                    signature_elem.set("desc", match_data["entry"]['description'])
                    signature_elem.set("name", match_data['entry']['short_name'])
                    if match_data['entry']['accession'] != "-":
                        entry_elem = ET.SubElement(signature_elem, "entry")
                        entry_elem.set("ac", match_data['entry']['accession'])
                        entry_elem.set("desc", match_data['entry']['description'])
                        entry_elem.set("name", match_data['entry']['short_name'])
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
                    if match_data['member_db'].upper() == "CDD":
                        location_elem = ET.SubElement(locations_elem, "cdd-location")
                        location_elem.set("end", str(location["end"]))
                        location_elem.set("start", str(location["start"]))
                        location_elem.set("representative", str(location["representative"]))
                        location_elem.set("evalue", str(location["evalue"]))
                        location_elem.set("score", str(location["score"]))
                        location_elem.set("postProcessed", str(location["postProcessed"]))
                    elif match_data['member_db'].upper() == "SIGNALP":
                        location_elem.set("end", str(location["end"]))
                        location_elem.set("start", str(location["start"]))
                        location_elem.set("pvalue", str(location["pvalue"]))
                    elif match_data['member_db'].upper() == "SFLD":
                        location_elem.set("sites", location["sites"])
                    else:
                        location_elem = ET.SubElement(locations_elem, "analysis-location")
                        location_elem.set("env-end", str(location["envelopeEnd"]))
                        location_elem.set("env-start", str(location["envelopeStart"]))
                        location_elem.set("score", str(location["score"]))
                        location_elem.set("evalue", str(location["evalue"]).upper())
                        location_elem.set("hmm-start", str(location["hmmStart"]))
                        location_elem.set("hmm-end", str(location["hmmEnd"]))
                        location_elem.set("hmm-length", str(location["hmmLength"]))
                        location_elem.set("hmm-bounds", str(location["hmmBounds"]))
                        location_elem.set("start", str(location["start"]))
                        location_elem.set("end", str(location["end"]))
                        location_elem.set("representative", str(location["representative"]))
                        if match_data['member_db'].upper() in ["NCBIFAM", "ANTIFAM"]:
                            location_elem.set("post-processed", str(location["postProcessed"]))
                        try:
                            location_elem.set("alignment", str(location["alignment"]))
                            location_elem.set("cigar-alignment", str(location["cigar_alignment"]))
                        except KeyError:
                            pass

                    location_frags_elem = ET.SubElement(location_elem, "location-fragments")
                    if 'sites' in location:
                        for site in location['sites']:
                            if match_data['member_db'].upper() == "CDD":
                                location_frag_elem = ET.SubElement(location_frags_elem, "cdd-location-fragment")
                                location_frag_elem.set("start", str(site['siteLocations']["start"]))
                                location_frag_elem.set("end", str(site['siteLocations']["end"]))
                                location_frag_elem.set("residue", str(site['siteLocations']["residue"]))
                            else:
                                for sitelocation in site['siteLocations']:
                                    location_frag_elem = ET.SubElement(location_frags_elem, "analysis-location-fragment")
                                    location_frag_elem.set("start", str(sitelocation["start"]))
                                    location_frag_elem.set("end", str(sitelocation["end"]))
                                    location_frag_elem.set("dc-status", "")

    tree = ET.ElementTree(root)
    ET.indent(tree, space="\t", level=0)
    tree.write(xml_output, encoding='utf-8')


def write_results(sequences_path: str, matches_path: str, output_format: list, output_path: str, version: str):
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
        tsv_output(seq_matches, output_path)
    if "TSV-PRO" in output_format:
        tsv_pro_output(seq_matches, output_path)
    if "JSON" in output_format:
        json_output(seq_matches, output_path, version)
    if "XML" in output_format:
        xml_output(seq_matches, output_path, version)


def main():
    args = sys.argv[1:]

    sequences = args[0]
    matches = args[1]
    formats_str = args[2]
    output_path = args[3]
    version = args[4]

    formats = formats_str.upper().split(',')
    write_results(sequences, matches, formats, output_path, version)


if __name__ == "__main__":
    main()
