import os
import xml.etree.ElementTree as ET


MATCH_ELEMENT = {
    'SIGNALP': 'signal-peptide',
    'CDD': 'cdd-domain',
    'ANTIFAM': 'hmmer3-match',
    'FUNFAM': 'hmmer3-match',
    'GENE3D': 'hmmer3-match',
    'HAMAP': 'hmmer3-match',
    'NCBIFAM': 'hmmer3-match',
    'PANTHER': 'hmmer3-match',
    'SFLD': 'hmmer3-match',
}

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
                        location_elem = ET.SubElement(locations_elem, "analysis-location")
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
                        location_elem = ET.SubElement(locations_elem, "analysis-location")
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
                                location_frag_elem = ET.SubElement(location_frags_elem, "analysis-location-fragment")
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
