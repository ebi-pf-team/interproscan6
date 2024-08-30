"""Test the python script that writes the XML output.

Fixtures are prefixed with 'x_' to indicate they are used in unit
testing the XML format writer.

These test are intened to be run from the root of the repository using:
python -m pytest -v
"""

import json
import xml.etree.ElementTree as ET

import pytest

from interproscan.scripts.output.format_writer import xml_output


@pytest.fixture
def x_expected_build_protein(test_output_dir):
    _path = test_output_dir / "format_writer/build_funcs/expected_xml_protein.xml"
    with open(_path, 'r') as f:
        tree = ET.parse(f)
    return xml_to_dict(tree.getroot())


@pytest.fixture
def x_expected_build_nucleic(test_output_dir):
    _path = test_output_dir / "format_writer/build_funcs/expected_xml_nucleic.xml"
    with open(_path, 'r') as f:
        tree = ET.parse(f)
    return xml_to_dict(tree.getroot())


@pytest.fixture
def x_matches_input_dir(test_input_dir):
    return test_input_dir / "format_writer/matches"


@pytest.fixture
def x_matches_output_dir(test_output_dir):
    return test_output_dir / "format_writer/matches/xml/"


@pytest.fixture
def x_nucleic_seq_matches(test_input_dir):
    _path = test_input_dir / "format_writer/nucleic/seq_match_dict.json"
    with open(_path, "r") as fh:
        _input = json.load(fh)
    return _input


@pytest.fixture
def x_prot_seq_matches(test_input_dir):
    _path = test_input_dir / "format_writer/protein/seq_match_dict.json"
    with open(_path, "r") as fh:
        _input = json.load(fh)
    return _input


@pytest.fixture
def x_protein_elm():
    root = ET.Element("protein-matches", xmlns="")
    root.set("interproscan-version", "unit.test")
    protein_elem = ET.SubElement(root, "protein")
    return protein_elem, root


@pytest.fixture
def x_outpath(test_output_dir):
    outdir = test_output_dir / "temp"
    outdir.mkdir(parents=True, exist_ok=True)
    filepath = outdir / "unittest.ips6.xml"
    return filepath


def load_match_data(member_db, matches_input_dir):
    with open((matches_input_dir / f"{member_db.upper()}.match-data.json"), "r") as fh:
        match_data = json.load(fh)
    return match_data


def load_expected_match_data(member_db, matches_output_dir):
    with open((matches_output_dir / f"{member_db.upper()}.matches.xml"), 'r') as f:
        tree = ET.parse(f)
    raw_dict = xml_to_dict(tree.getroot())
    return parse_expected_matches_dict(raw_dict['protein'][0]['matches'][0])


def parse_expected_matches_dict(matches: list[dict[str, dict]]) -> dict[str, dict[str, dict]]:
    match_dict = {}
    for domain_type, signature_list in matches.items():
        for signature in signature_list:
            sig_acc = signature['signature'][0]['attributes']['ac']
            match_dict[sig_acc] = {
                'signature': {
                    'acc': sig_acc,
                    'name': signature['signature'][0]['attributes']['name'],
                    'description': signature['signature'][0]['attributes']['desc'],
                },
                'signature-library-release': {
                    'library': signature['signature'][
                        0]['signature-library-release'][0]['attributes']['library'],
                    'version': signature['signature'][
                        0]['signature-library-release'][0]['attributes']['version'],
                },
                'model-ac': signature['model-ac'][0],
                'locations': []
            }
            for location in signature['locations'][0]:
                new_location = {}
                for key in location['analysis-location'][0]['attributes']:
                    new_location[key] = location['analysis-location'][0]['attributes'][key]
                match_dict[sig_acc]['locations'].append(new_location)

    return match_dict


def xml_to_dict(element: ET.Element):
    result = {}
    if element.attrib:
        result.update({'attributes': element.attrib})
    for child in element:
        if child.tag not in result:
            result[child.tag] = []
        result[child.tag].append(xml_to_dict(child))
    for key, value in result.items():
        if isinstance(value, list) and all(isinstance(item, dict) for item in value):
            result[key] = sorted(value, key=lambda x: str(x))
    return result


def test_build_xml_protein(x_prot_seq_matches, x_outpath, x_expected_build_protein, monkeypatch):
    def mock_get_matches(*args, **kwards):
        root = ET.Element("protein-matches", xmlns="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6/schemas")
        root.set("interproscan-version", "unit.test")
        protein_elem = ET.SubElement(root, "protein")
        matches_elem = ET.SubElement(protein_elem, "matches")
        return protein_elem

    monkeypatch.setattr(xml_output, "add_xml_output_matches", mock_get_matches)

    xml_output.build_xml_output_protein(x_prot_seq_matches, x_outpath, "unit.test")
    with open(x_outpath, 'r') as f:
        result_tree = ET.parse(f)
    result_tree = result_tree.getroot()

    assert xml_to_dict(result_tree) == x_expected_build_protein, "Mistmatch in overall XML when input is protein sequences"

    x_outpath.unlink()


def test_build_xml_nucleic(x_nucleic_seq_matches, x_outpath, x_expected_build_nucleic, monkeypatch):
    def mock_get_matches(*args, **kwards):
        root = ET.Element("nucleic-matches", xmlns="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6/schemas")
        root.set("interproscan-version", "unit.test")
        nt_seq_elem = ET.SubElement(root, "nucleotide-sequence")
        sequence_elem = ET.SubElement(nt_seq_elem, "sequence")
        xref_elem = ET.SubElement(nt_seq_elem, "xref")
        orf_elem = ET.SubElement(nt_seq_elem, "orf")
        protein_elem = ET.SubElement(orf_elem, "protein")
        sequence_elem = ET.SubElement(protein_elem, "sequence")
        xref_elem = ET.SubElement(protein_elem, "xref")
        matches_elem = ET.SubElement(protein_elem, "matches")
        return protein_elem

    monkeypatch.setattr(xml_output, "add_xml_output_matches", mock_get_matches)

    xml_output.build_xml_output_nucleic(x_nucleic_seq_matches, x_outpath, "unit.test")
    with open(x_outpath, 'r') as f:
        result_tree = ET.parse(f)
    result_tree = result_tree.getroot()

    assert xml_to_dict(result_tree) == x_expected_build_nucleic, "Mistmatch in overall XML when input is protein sequences"

    x_outpath.unlink()


def test_cdd_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "CDD"
    match_data = load_match_data(member_db, x_matches_input_dir)

    result = xml_output.add_xml_output_matches(x_protein_elm[0], match_data)
    result_dict = xml_to_dict(result)
    expected_dict = load_expected_match_data(member_db, x_matches_output_dir)
    


    # tree = ET.ElementTree(x_protein_elm[1])
    # ET.indent(tree, space="\t", level=0)
    # tree.write((x_matches_output_dir / f"{member_db}.matches.xml"), encoding='utf-8')
