import xml.etree.ElementTree as ET

from .regex import NT_KEY_PATTERN, NT_SEQ_ID_PATTERN


MATCH_ELEMENT = {
    'ANTIFAM': 'hmmer3-match',
    'CDD': 'cdd-domain',
    'COILS': 'coils',
    'FUNFAM': 'hmmer3-match',
    'GENE3D': 'hmmer3-match',
    'HAMAP': 'hmmer3-match',
    'MOBIDB': 'mobidb-match',
    'NCBIFAM': 'hmmer3-match',
    'PANTHER': 'hmmer3-match',
    'PFAM': 'hmmer3-match',
    'PHOBIUS': 'phobius-match',
    'PIRSF': 'hmmer3-match',
    'PIRSR': 'hmmer3-match',
    'PRINTS': 'fingerprints-match',
    'PROSITE_PATTERNS': 'profilescan-match',
    'PROSITE_PROFILES': 'profilesearch-match',
    'SFLD': 'hmmer3-match',
    'SIGNALP': 'signalp-match',
    'SIGNALP_EUK': 'signalp-euk-match',
    'SMART': 'hmmer2-match',
    'SUPERFAMILY': 'hmmer3-match',
}


def build_xml_output_protein(seq_matches: dict, output_path: str, version: str):
    """Build the root of the XML when the input to IPS6 is Protein sequences"""
    root = ET.Element("protein-matches", xmlns="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6/schemas")
    root.set("interproscan-version", version)

    for seq_id, data in seq_matches.items():
        protein_elem = ET.SubElement(root, "protein")
        sequence_elem = ET.SubElement(protein_elem, "sequence")
        sequence_elem.text = data['sequences']['sequence']
        sequence_elem.set("md5", data['sequences']['md5'])
        xref_elem = ET.SubElement(protein_elem, "xref")
        xref_elem.set("id", seq_id)
        xref_elem.set("name", data['sequences']['seq_id'])
        protein_elem = add_xml_output_matches(protein_elem, data)

    tree = ET.ElementTree(root)
    ET.indent(tree, space="\t", level=0)
    tree.write(output_path, encoding='utf-8')


def build_xml_output_nucleic(seq_matches: dict, output_path: str, version: str):
    """Build the root of the XML when the input to IPS6 is nucleic sequences"""
    def _get_orf_keys(seq_matches: dict):
        """Retrieve the keys from seq_matches for all ORFs for each input nucleic seq"""
        current_nt = ""
        nt_keys = {}
        for key in seq_matches:
            nt = NT_KEY_PATTERN.match(key).group(1)
            if nt != current_nt:
                current_nt = nt
                nt_keys[nt] = [key]
            else:
                nt_keys[nt].append(key)
        return nt_keys

    root = ET.Element("nucleotide-sequence-matches", xmlns="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6/schemas")
    root.set("interproscan-version", version)

    nt_keys = _get_orf_keys(seq_matches)
    for nt_seq_id, orf_ids in nt_keys.items():
        # grab the data for the nucleic seq from the first ORF
        data = seq_matches[orf_ids[0]]
        nt_seq_elem = ET.SubElement(root, "nucleotide-sequence")
        sequence_elem = ET.SubElement(nt_seq_elem, "sequence")
        sequence_elem.text = data['sequences']['nt_sequence']
        sequence_elem.set("md5", data['sequences']['nt_md5'])
        xref_elem = ET.SubElement(nt_seq_elem, "xref")
        xref_elem.set("id", data['sequences']['nt_seq_id'].split(maxsplit=1)[0])
        xref_elem.set("name", data['sequences']['nt_seq_id'])

        # add data for each ORF. One ORF = One Protein elem
        for seq_id in orf_ids:
            data = seq_matches[seq_id]
            orf_elem = ET.SubElement(nt_seq_elem, "orf")
            seq_data = NT_SEQ_ID_PATTERN.match(data['sequences']['seq_id'])
            strand = "SENSE" if int(seq_data.group(3)) < 4 else "ANTISENSE"
            coords = seq_data.group(2).split("..")
            orf_elem.set("end", coords[1])
            orf_elem.set("start", coords[0])
            orf_elem.set("strand", strand)

            protein_elem = ET.SubElement(orf_elem, "protein")
            sequence_elem = ET.SubElement(protein_elem, "sequence")
            sequence_elem.text = data['sequences']['sequence']
            sequence_elem.set("md5", data['sequences']['md5'])
            xref_elem = ET.SubElement(protein_elem, "xref")
            xref_elem.set("id", f"orf{seq_id.split('_orf')[1]}")
            xref_elem.set("name", data['sequences']['seq_id'])
            protein_elem = add_xml_output_matches(protein_elem, data)

    tree = ET.ElementTree(root)
    ET.indent(tree, space="\t", level=0)
    tree.write(output_path, encoding='utf-8')


def add_xml_output_matches(protein_elem: ET.SubElement, data: dict):
    """Add the matches to the XML tree"""
    def _check_null(value, acc=False):
        return str(value) if str(value).lower() not in ["none", "null", ""] else "-" if acc else ""

    matches_elem = ET.SubElement(protein_elem, "matches")
    if 'matches' in data and data['matches']:
        for match_key, match_data in data['matches'].items():

            # Write each signal peptide/transM domain/(non-)cytoplasmic location as
            # a separate phobius-match
            if match_data['member_db'].upper() == "PHOBIUS":
                for location in match_data["locations"]:
                    match_elem = ET.SubElement(matches_elem, MATCH_ELEMENT[match_data['member_db'].upper()])
                    signature_elem = ET.SubElement(match_elem, "signature")
                    signature_elem.set("ac", match_data['accession'])
                    signature_elem.set("desc", _check_null(match_data['entry']['description']))
                    signature_elem.set("name", _check_null(match_data['entry']['name']))
                    signature_library_elem = ET.SubElement(signature_elem, "signature-library-release")
                    signature_library_elem.set("library", match_data['member_db'].upper())
                    signature_library_elem.set("version", match_data['version'])
                    model_ac_elem = ET.SubElement(match_elem, "model-ac")
                    model_ac_elem.text = match_key
                    locations_elem = ET.SubElement(match_elem, "locations")
                    location_elem = ET.SubElement(locations_elem, "analysis-location")
                    location_elem.set("start", str(location["start"]))
                    location_elem.set("end", str(location["end"]))
                    location_elem.set("representative", str(location["representative"]))
                    location_frags_elem = ET.SubElement(location_elem, "location-fragments")
                    location_frag_elem = ET.SubElement(location_frags_elem, "analysis-location-fragment")
                    location_frag_elem.set("start", str(location["location-fragments"][0]["start"]))
                    location_frag_elem.set("end", str(location["location-fragments"][0]["end"]))
                    location_frag_elem.set("dc-status", str(location["location-fragments"][0]["dc-status"]))
                continue

            match_elem = ET.SubElement(matches_elem, MATCH_ELEMENT[match_data['member_db'].upper()])

            # prevent evalue and score from MLS appearing in xml format for CDD matches
            if match_data['member_db'].upper() != "CDD":
                try:
                    match_elem.set("evalue", str(match_data['evalue']).upper())
                    match_elem.set("score", str(match_data["score"]))
                except KeyError:
                    pass  # some members may not have evalue or score on this level (e.g. cdd)

            if match_data['member_db'].upper() == "PANTHER":
                match_elem.set("protein-class", _check_null(match_data['proteinClass']))
                match_elem.set("graft-point", _check_null(match_data['graftPoint']))

            signature_elem = ET.SubElement(match_elem, "signature")
            if match_data['member_db'].upper() not in ['SIGNALP', 'SIGNALP_EUK', 'SUPERFAMILY']:
                # member db that don't have sigs, so no accs etc.
                signature_elem.set("ac", match_data['accession'])

            if match_data['member_db'].upper() in ['SIGNALP', 'SIGNALP_EUK']:
                signature_elem.set("ac", match_data['accession'])
                signature_elem.set("name", match_data['name'])

            if match_data['entry']:
                acc = _check_null(match_data['entry']['accession'], acc=True)
                sig_type = _check_null(match_data['entry']['type'])
                desc = _check_null(match_data["entry"]['description'])
                name = _check_null(match_data["entry"]['name'])
                entry_desc = _check_null(match_data["entry"].get('ipr_description', "-"))
                entry_name = _check_null(match_data["entry"].get('ipr_name', "-"))
                signature_elem.set("desc", desc)
                signature_elem.set("name", name)
                if match_data['entry']['accession'] != "-":
                    entry_elem = ET.SubElement(signature_elem, "entry")
                    entry_elem.set("ac", acc)
                    entry_elem.set("desc", entry_desc)
                    entry_elem.set("name", entry_name)
                    entry_elem.set("type", sig_type)
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
            if match_data['entry']:
                signature_library_elem.set("version", match_data['version'])
            else:
                signature_library_elem.set("version", "-")
            model_ac_elem = ET.SubElement(match_elem, "model-ac")
            model_ac_elem.text = match_key

            locations_elem = ET.SubElement(match_elem, "locations")

            for location in match_data['locations']:
                if match_data['member_db'].upper() == "CDD":
                    location_elem = ET.SubElement(locations_elem, "analysis-location")
                    location_elem.set("end", str(location["end"]))
                    location_elem.set("start", str(location["start"]))
                    location_elem.set("representative", str(location["representative"]))
                    location_elem.set("evalue", str(location["evalue"]))
                    location_elem.set("score", str(location["score"]))

                elif match_data['member_db'].upper() == "HAMAP":
                    location_elem = ET.SubElement(locations_elem,"analysis-location")
                    location_elem.set("score", str(location["score"]))
                    location_elem.set("alignment", str(location["alignment"]))

                elif match_data['member_db'].upper() == "PRINTS":
                    location_elem = ET.SubElement(locations_elem, "analysis-location")
                    location_elem.set("motifNumber", str(int(location["motifNumber"])))
                    location_elem.set("pvalue", str(location["pvalue"]))
                    location_elem.set("score", str(location["score"]))
                    location_elem.set("end", str(location["end"]))
                    location_elem.set("start", str(location["start"]))
                    location_elem.set("representative", str(location["representative"]))

                elif match_data['member_db'].upper() == "PROSITE_PROFILES":
                    location_elem = ET.SubElement(locations_elem, "analysis-location")
                    location_elem.set("score", str(location["score"]))
                    location_elem.set("start", str(location["start"]))
                    location_elem.set("end", str(location["end"]))
                    location_elem.set("alignment", str(location["alignment"]))

                elif match_data['member_db'].upper() == "PROSITE_PATTERNS":
                    location_elem = ET.SubElement(locations_elem, "analysis-location")
                    location_elem.set("start", str(location["start"]))
                    location_elem.set("end", str(location["end"]))
                    location_elem.set("alignment", str(location["alignment"]))
                    location_elem.set("cigar-alignment", str(location["cigarAlignment"]))

                elif match_data['member_db'].upper() in ["SIGNALP", "SIGNALP_EUK"]:
                    location_elem = ET.SubElement(locations_elem,"analysis-location")
                    location_elem.set("end", str(location["end"]))
                    location_elem.set("start", str(location["start"]))
                    location_elem.set("pvalue", str(location["pvalue"]))
                    location_elem.set("cleavage_start", str(location["cleavage_start"]))
                    location_elem.set("cleavage_end", str(location["cleavage_end"]))
                    location_elem.set("representative", str(location["representative"]))

                elif match_data['member_db'].upper() == "SFLD":
                    location_elem = ET.SubElement(locations_elem, "analysis-location")
                    location_elem.set("evalue", str(location["evalue"]).upper())
                    location_elem.set("score", str(location["score"]))
                    location_elem.set("hmm-start", str(location["hmmStart"]))
                    location_elem.set("hmm-end", str(location["hmmEnd"]))
                    location_elem.set("env-end", str(location["envelopeEnd"]))
                    location_elem.set("env-start", str(location["envelopeStart"]))
                    try:
                        location_elem.set("sites", location["sites"])
                    except KeyError:
                        location_elem.set("sites", [])

                elif match_data['member_db'].upper() == "PIRSR":
                    location_elem = ET.SubElement(locations_elem, "analysis-location")
                    try:
                        location_elem.set("sites", location["sites"])
                    except KeyError:
                        location_elem.set("sites", [])

                elif match_data['member_db'].upper() == "SMART":
                    location_elem = ET.SubElement(locations_elem, "analysis-location")
                    location_elem.set("score", str(location["score"]))
                    location_elem.set("evalue", str(location["evalue"]))
                    location_elem.set("hmm-start", str(location["hmmStart"]))
                    location_elem.set("hmm-end", str(location["hmmEnd"]))
                    location_elem.set("hmm-len", str(location["hmmLength"]))
                    location_elem.set("hmm-bounds", str(location["hmmBounds"]))
                    location_elem.set("start", str(location["start"]))
                    location_elem.set("end", str(location["end"]))
                    location_elem.set("representative", str(location["representative"]))

                elif match_data['member_db'].upper() in ["COILS", "MOBIDB", "SUPERFAMILY"]:
                    location_elem = ET.SubElement(locations_elem, "analysis-location")
                    location_elem.set("start", str(location["start"]))
                    location_elem.set("end", str(location["end"]))
                    location_elem.set("representative", str(location["representative"]))

                else:
                    location_elem = ET.SubElement(locations_elem, "analysis-location")
                    location_elem.set("env-end", str(location["envelopeEnd"]))
                    location_elem.set("env-start", str(location["envelopeStart"]))
                    location_elem.set("score", str(location["score"]))
                    location_elem.set("evalue", str(location["evalue"]).upper())

                    if match_data['member_db'].upper() == 'PIRSF':
                        # PIRSF uses the ali from/to for the hmmstart/end
                        location_elem.set("hmm-start", str(location["start"]))
                        location_elem.set("hmm-end", str(location["end"]))
                    else:
                        location_elem.set("hmm-start", str(location["hmmStart"]))
                        location_elem.set("hmm-end", str(location["hmmEnd"]))

                    location_elem.set("hmm-length", str(location["hmmLength"]))
                    location_elem.set("hmm-bounds", str(location["hmmBounds"]))

                    if match_data['member_db'].upper() == 'PIRSF':
                        # PIRSF uses the env from/to for the start/end
                        location_elem.set("start", str(location["envelopeStart"]))
                        location_elem.set("end", str(location["envelopeEnd"]))
                    else:
                        location_elem.set("start", str(location["start"]))
                        location_elem.set("end", str(location["end"]))

                    location_elem.set("representative", str(location["representative"]))
                    try:
                        location_elem.set("alignment", str(location["alignment"]))
                        location_elem.set("cigar-alignment", str(location["cigar_alignment"]))
                    except KeyError:
                        pass

                location_frags_elem = ET.SubElement(location_elem, "location-fragments")
                if 'sites' in location:
                    for site in location['sites']:
                        if match_data['member_db'].upper() == "CDD":
                            location_frag_elem = ET.SubElement(
                                location_frags_elem,
                                "analysis-location-fragment"
                            )
                            location_frag_elem.set("description", str(site['description']))
                            location_frag_elem.set("numLocations", str(site['numLocations']))
                        else:
                            for sitelocation in site['siteLocations']:
                                location_frag_elem = ET.SubElement(
                                    location_frags_elem,
                                    "analysis-location-fragment"
                                )
                                location_frag_elem.set("start", str(sitelocation["start"]))
                                location_frag_elem.set("end", str(sitelocation["end"]))
                                location_frag_elem.set("residue", str(sitelocation["residue"]))
                if 'location-fragments' in location:
                    for location_fragment in location['location-fragments']:
                        location_frag_elem = ET.SubElement(location_frags_elem, "analysis-location-fragment")
                        location_frag_elem.set("start", str(location_fragment["start"]))
                        location_frag_elem.set("end", str(location_fragment["end"]))
                        location_frag_elem.set("dc-status", str(location_fragment["dc-status"]))

    return protein_elem