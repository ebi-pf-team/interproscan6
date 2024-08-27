import json
import os

from .regex import NT_SEQ_ID_PATTERN


BOOLEAN_MAP = {"true": True, "false": False}


def build_json_output_protein(seq_matches: dict, output_path: str, version: str):
    results = []

    for seq_id, data in seq_matches.items():
        results.append({
            "sequence": data['sequences']['sequence'],
            "md5": data['sequences']['md5'],
            "matches": get_matches(data),
            "xref": [{
                "name": data['sequences']['seq_id'],
                "id": seq_id
            }]
        })

    final_data = {"interproscan-version": version, 'results': results}
    with open(output_path, 'w') as json_file:
        json_file.write(json.dumps(final_data, indent=2))

    return final_data


def build_json_output_nucleic(seq_matches: dict, output_path: str, version: str):
    """Iterate through the ORFs in seq_matches."""
    json_output = os.path.join(output_path + '.json')
    results = {}  # nucleic seq id: {}

    # seq_id = <nucleic acid seq id>_orf<id>
    for seq_id, data in seq_matches.items():  # iterates the ORFs
        nucleic_seq_id = data['sequences']['nt_seq_id'].split(maxsplit=1)[0]

        # add nucleic seq to the results
        if nucleic_seq_id not in results:
            results[nucleic_seq_id] = {
                "sequence": data['sequences']['nt_sequence'],
                "md5": data['sequences']['nt_md5'],
                "crossReferences": [{
                    'name': data['sequences']['nt_seq_id'],
                    'id': nucleic_seq_id,
                }],
                "openReadingFrames": []
            }

        # create and add a dict for the open reading frame
        nt_match = NT_SEQ_ID_PATTERN.match(data['sequences']['seq_id'])
        results[nucleic_seq_id]["openReadingFrames"].append(
            {
                "start": int(nt_match.group(2).split("..")[0]),
                "end": int(nt_match.group(2).split("..")[1]),
                "strand": "SENSE" if int(nt_match.group(3)) < 4 else "ANTISENSE",
                "protein": {
                    "sequence": data['sequences']['sequence'],
                    "md5": data['sequences']['md5'],
                    "matches": get_matches(data),
                    "xref": [{
                        "name": data['sequences']['seq_id'],
                        "id": data['sequences']['seq_id'].split(maxsplit=1)[0],
                    }]
                },
            }
        )

    final_data = {"interproscan-version": version, 'results': list(results.values())}
    with open(output_path, 'w') as json_file:
        json_file.write(json.dumps(final_data, indent=2))

    return final_data


def get_matches(data: dict):
    matches = []

    if 'matches' in data and data['matches']:
        for match_key, match_data in data['matches'].items():  # match_key == sig_Acc
            try:
                description = match_data['description']
            except KeyError:
                description = None if match_data['member_db'].upper() == "SUPERFAMILY" else "-"

            entry = None
            if match_data['entry']:
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

            # drop the subfamily, else if gene3d/funfam then keep info after ":", e.g. G3DSA:3.20.20.70
            accession = match_data['accession'].split(":")[0] \
                if match_data['member_db'].upper() not in ["GENE3D", "FUNFAM"] \
                else match_data['accession']

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

            # For these member dbs we write each domain location as a separate
            # 'signature' match in the final results
            if match_data['member_db'].upper() in ["MOBIDB", "PHOBIUS", "SUPERFAMILY"]:
                for location in match_data['locations']:
                    location_info = {
                        "start": int(location["start"]),
                        "end": int(location["end"]),
                        "representative": BOOLEAN_MAP.get(location["representative"].lower(), False),
                        "location-fragments": [{
                            "start": int(location["start"]),
                            "end": int(location["end"]),
                            "dc-status": "CONTINUOUS"
                        }]
                    }
                    if match_data['member_db'].upper() == "PHOBIUS":
                        match = {
                            "signature": signature,
                            "locations": [location_info],
                            "model-ac": accession
                        }

                    elif match_data['member_db'].upper() == "SUPERFAMILY":
                        location_info['evalue'] = float(location['evalue'])
                        try:
                            location_info["hmmLength"] = match_data['hmm_length']
                        except KeyError:
                            location_info["hmmLength"] = location['hmmLength']
                        match = {
                            "signature": signature,
                            "locations": [location_info],
                            "evalue": float(match_data["evalue"]),
                            "model-ac": match_data.get('model-ac', match_data['accession'])
                        }
                    else:
                        location_info["sequence-feature"] = location["sequence-feature"]
                        match = {
                            "signature": signature,
                            "locations": [location_info]
                        }

                    matches.append(match)
            else:
                if len(match_data['locations']) > 0:
                    locations = []
                    for location in match_data['locations']:
                        # PIRSF usse the envelope start and stop
                        if match_data['member_db'].upper() == "PIRSF":
                            info = {
                                "start": int(location["envelopeStart"]),
                                "end": int(location["envelopeEnd"])
                            }
                        else:
                            info = {
                                "start": int(location["start"]),
                                "end": int(location["end"])
                            }

                        info["representative"] = BOOLEAN_MAP.get(
                            location["representative"].lower(),
                            False
                        )

                        if match_data['member_db'].upper() == "CDD":
                            info["evalue"] = float(location["evalue"])
                            info["score"] = float(location["score"])

                        elif match_data['member_db'].upper() in ["COILS", "PHOBIUS"]:
                            pass  # data alreadt listed in into

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

                        elif match_data['member_db'].upper() == "PIRSF":
                            # PIRSF uses the ali from (start) and ali to (end)
                            # for the hmmStart and hmmEnd
                            # and env from/to for the start and end
                            info["evalue"] = float(location["evalue"])
                            info["score"] = float(location["score"])
                            info["hmmStart"] = int(location["start"])
                            info["hmmEnd"] = int(location["end"])
                            info["hmmLength"] = int(location["hmmLength"])
                            info["hmmBounds"] = location["hmmBounds"]
                            info["envelopeStart"] = int(location["envelopeStart"])
                            info["envelopeEnd"] = int(location["envelopeEnd"])

                        elif match_data['member_db'].upper() == "PRINTS":
                            info["pvalue"] = float(location["pvalue"])
                            info["score"] = float(location["score"])
                            info["motifNumber"] = int(location["motifNumber"])

                        elif match_data['member_db'].upper() == "PROSITE_PROFILES":
                            info["score"] = float(location["score"])
                            info["alignment"] = str(location["alignment"])

                        elif match_data['member_db'].upper() == "PROSITE_PATTERNS":
                            info["cigarAlignment"] = location["cigarAlignment"]
                            info["alignment"] = location["alignment"]
                            info["level"] = location["level"]

                        elif match_data['member_db'].upper() in ["PIRSR", "SFLD"]:
                            info["evalue"] = float(location["evalue"])
                            info["score"] = float(location["score"])
                            info["hmmStart"] = int(location["hmmStart"])
                            info["hmmEnd"] = int(location["hmmEnd"])
                            info["hmmLength"] = int(location["hmmLength"])
                            info["envelopeStart"] = int(location["envelopeStart"])
                            info["envelopeEnd"] = int(location["envelopeEnd"])

                        elif match_data['member_db'].upper() in ["SIGNALP", "SIGNALP_EUK"]:
                            info["pvalue"] = float(location["pvalue"])
                            info["cleavageStart"] = location["cleavage_start"]
                            info["cleavageEnd"] = location["cleavage_end"]

                        elif match_data['member_db'].upper() == "SMART":
                            info["evalue"] = float(location["evalue"])
                            info["score"] = float(location["score"])
                            info["hmmStart"] = int(location["hmmStart"])
                            info["hmmEnd"] = int(location["hmmEnd"])
                            info["hmmLength"] = int(location["hmmLength"])
                            info["hmmBounds"] = location["hmmBounds"]

                        else:
                            info["evalue"] = float(location["evalue"])
                            info["score"] = float(location["score"])
                            info["hmmStart"] = int(location["hmmStart"])
                            info["hmmEnd"] = int(location["hmmEnd"])
                            info["hmmLength"] = int(location["hmmLength"])
                            info["hmmBounds"] = location["hmmBounds"]
                            info["envelopeStart"] = int(location["envelopeStart"])
                            info["envelopeEnd"] = int(location["envelopeEnd"])

                        if match_data['member_db'].upper() in ["CDD", "PIRSR", "SFLD"]:
                            info["sites"] = location["sites"] if "sites" in location else []

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

                    if match_data['member_db'].upper() not in [
                        "CDD", "COILS", "HAMAP", "PHOBIUS", "PIRSR",
                        "PROSITE_PROFILES", "PROSITE_PATTERNS",
                        "PRINTS", "SIGNALP", "SIGNALP_EUK"
                    ]:
                        match["evalue"] = float(match_data['evalue'])
                        match["score"] = float(match_data['score'])

                    if 'model-ac' in match_data:
                        match["model-ac"] = match_data['model-ac']
                    else:
                        match["model-ac"] = match_data['accession']

                    if match_data['member_db'].upper() == "SFLD":
                        match["scope"] = None

                    elif match_data['member_db'].upper() == "PANTHER":
                        match["name"] = match_data['entry']['family_name']
                        match['accession'] = match_data['accession']
                        match["goXRefs"] = entry["goXRefs"] if entry else []
                        signature["description"] = None
                        signature["name"] = match_data['entry']['description']
                        match['proteinClass'] = match_data['proteinClass']
                        match['graftPoint'] = match_data['graftPoint']

                    elif match_data['member_db'].upper() == "PRINTS":
                        match["evalue"] = float(match_data['evalue'])
                        match["graphscan"] = str(match_data["graphscan"])

                    elif match_data['member_db'].upper() in ["SIGNALP", "SIGNALP_EUK"]:
                        match["orgType"] = match_data["orgType"]

                    matches.append(match)

    return matches
