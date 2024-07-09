import json
import os


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

                elif match_key == "transmembrane_prediction":
                    locations = []
                    for location in match_data['locations']:
                        info = {"location_tag": location["location_tag"],
                            "start": int(location["start"]),
                            "end": int(location["end"])
                        }
                        locations.append(info)
                    match = {"signature": match_key,
                        "DeepTMHMM_release": match_data["version"],
                        "locations": locations
                    }

                else:
                    try:
                        description = match_data['description']
                    except KeyError:
                        description = "-"
                    entry = None
                    if match_data['entry']['accession']:
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
                            pass
                        try:
                            entry["pathwayXRefs"] = match_data['entry']['pathwayXRefs']
                        except KeyError:
                            pass

                    if match_data['member_db'].upper() == "CDD":
                        description = match_data['name']

                    if match_data['member_db'].upper() in ["GENE3D", "FUNFAM"]:
                        accession = match_data['accession']  # GENE3D needs the info after ":" (e.g G3DSA:3.20.20.70)
                    else:
                        accession = match_data['accession'].split(":")[0]  # drop subfamily

                    signature = {
                        "accession": accession,
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
                            "representative": boolean_map.get(location["representative"].lower(), False)
                        }
                        if match_data['member_db'].upper() == "CDD":
                            info["evalue"] = float(location["evalue"])
                            info["score"] = float(location["score"])

                        elif match_data['member_db'].upper() == "HAMAP":
                            info["score"] = float(location["score"])
                            info["alignment"] = location["alignment"]

                        elif match_data['member_db'].upper() == "PANTHER":
                            info["hmmStart"] = int(location["hmmStart"])
                            info["hmmEnd"] = int(location["hmmEnd"])
                            info["hmmLength"] = 0  # we have hmmLength but in i5 result its always 0
                            info["hmmBounds"] = location["hmmBounds"]
                            info["envelopeStart"] = int(location["envelopeStart"])
                            info["envelopeEnd"] = int(location["envelopeEnd"])

                        elif match_data['member_db'].upper() == "SFLD":
                            info["evalue"] = float(location["evalue"])
                            info["score"] = float(location["score"])
                            info["hmmStart"] = int(location["hmmStart"])
                            info["hmmEnd"] = int(location["hmmEnd"])
                            info["hmmLength"] = int(location["hmmLength"])
                            info["envelopeStart"] = int(location["envelopeStart"])
                            info["envelopeEnd"] = int(location["envelopeEnd"])

                        elif match_data['member_db'].upper() == "PROSITE_PROFILES":
                            info["score"] = float(location["score"])
                            info["alignment"] = str(location["alignment"])

                        else:
                            info["evalue"] = float(location["evalue"])
                            info["score"] = float(location["score"])
                            info["hmmStart"] = int(location["hmmStart"])
                            info["hmmEnd"] = int(location["hmmEnd"])
                            info["hmmLength"] = int(location["hmmLength"])
                            info["hmmBounds"] = location["hmmBounds"]
                            info["envelopeStart"] = int(location["envelopeStart"])
                            info["envelopeEnd"] = int(location["envelopeEnd"])
                            info["postProcessed"] = boolean_map.get(location["postProcessed"].lower())

                        if match_data['member_db'].upper() in ["SFLD", "CDD"]:
                            try:
                                info["sites"] = location["sites"]
                            except KeyError:
                                info["sites"] = []
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

                    if match_data['member_db'].upper() not in ["CDD", "HAMAP", "PROSITE_PROFILES"]:
                        match["evalue"] = float(match_data['evalue'])
                        match["score"] = float(match_data['score'])

                    if 'model-ac' in match_data:
                        match["model-ac"] = match_data['model-ac']
                    else:
                        match["model-ac"] = match_data['accession']

                    if match_data['member_db'].upper() == "SFLD":
                        match["scope"] = None

                    if match_data['member_db'].upper() == "PANTHER":
                        match["name"] = match_data['entry']['family_name']
                        match['accession'] = match_data['accession']
                        match["goXRefs"] = entry["goXRefs"] if entry else []
                        signature["description"] = None
                        signature["name"] = match_data['entry']['description']

                        # get protein class and graftpoint for Panther
                        try:
                            match['proteinClass'] = match_data['proteinClass']
                            match['graftPoint'] = match_data['graftPoint']
                        except KeyError:
                            pass

                if len(match_data['locations']) > 0:  # skip matches with no locations (we need to make sure it's valid to all members)
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
