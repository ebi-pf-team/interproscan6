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
    return parse_matches_dict(raw_dict['protein'][0]['matches'][0])


def parse_matches_dict(matches: list[dict[str, dict]]) -> dict[str, dict[str, dict]]:
    match_dict = {}
    for domain_type, signature_list in matches.items():
        for signature in signature_list:
            sig_acc = signature['signature'][0]['attributes']['ac'] if \
                'ac' in signature['signature'][0]['attributes'] else ''
            match_dict[sig_acc] = {
                'signature': {
                    'acc': sig_acc,
                    'name': signature['signature'][0]['attributes']['name'],
                    'description': \
                        signature['signature'][0]['attributes']['desc'] if 'desc' in \
                            signature['signature'][0]['attributes'] else '',
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
            for location in signature['locations']:
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

    assert xml_to_dict(result_tree) == x_expected_build_protein, \
        "Mistmatch in overall XML when input is protein sequences"

    # x_outpath.unlink()


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


def compare_signature_details(
    match_data: dict, expected_data: dict,
    sig_acc: str, member_db: str
) -> None:
    assert match_data['model-ac'] == expected_data['model-ac'], \
        f"Mismatches model-ac for {sig_acc}, {member_db}"
    assert all(
        match_data['signature'][key] == expected_data['signature'][key]
        for key in ['acc', 'name', 'description']
    ), f"Mismatch in 'signature' details for {sig_acc}, {member_db}"
    assert all(
        match_data['signature-library-release'][key] == expected_data['signature-library-release'][key]
        for key in ['library', 'version']
    ), f"Mismatch in 'signatureLibraryRelease' details for {sig_acc}, {member_db}"


def compare_location_details(
    locations: list[dict],
    expected_locations: list[dict],
    location_keys: list[str],
    sig_acc: str,
    member_db: str
) -> None:
    """
    :param locations: list == match_data['locations']
    :param expected_locations: list == expected_data['locations']
    :param location_keys: list, keys in 'locations' to compare
    :param sig_acc: str, accession of the associated signature
    :param member_db: str, name of the associated member db
    """
    assert len(locations) == len(expected_locations), f"Mismatched number of locations for {sig_acc}, {member_db}"
    for location in locations:
        found = False
        for expected_location in expected_locations:
            if location['start'] == expected_location['start'] \
                and location['end'] == expected_location['end']:
                found = True
                assert all(
                    location[key] == expected_location[key]
                    for key in location_keys
                ), (
                    f"Mismatch in 'location' details for {sig_acc}, {member_db} "
                    f"(start: {location['start']} - end {location['end']})"
                )
                if 'sites' in location:
                    assert len(location['sites']) == len(expected_location['sites'])
                if 'location-fragments' in location:
                    assert len(location['location-fragments']) == len(expected_location['location-fragments'])
                # gritty detailed comparison can be handeled by the intergration test
        if not found:
            pytest.fail((
                f"Location with start {location['start']} and end {location['end']} "
                f"not found in expected data for {sig_acc}, {member_db}"
            ))
    for expected_location in expected_locations:
        found = False
        for location in locations:
            if location['start'] == expected_location['start'] \
                and location['end'] == expected_location['end']:
                found = True
        if not found:
            pytest.fail((
                f"Expected location with start {location['start']} and end {location['end']} "
                f"not found in test result data for {sig_acc}, {member_db}"
            ))


# def test_antifam_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
#     member_db = "ANTIFAM"
#     match_data = load_match_data(member_db, x_matches_input_dir)

#     result = xml_output.add_xml_output_matches(x_protein_elm[0], match_data)
#     # # result_dict = xml_to_dict(result)
#     # result_dict = parse_matches_dict(result_dict['matches'][0])
#     # expected_dict = load_expected_match_data(member_db, x_matches_output_dir)

#     # for sig_acc, match_data in result_dict.items():
#     #     assert sig_acc in expected_dict, f"Signature {sig_acc} not in expected results, {member_db}"
#     #     if sig_acc in expected_dict:
#     #         expected_data = expected_dict[sig_acc]
#     #         compare_signature_details(match_data, expected_data, sig_acc, member_db)
#     #         locations = match_data['locations']
#     #         expected_locations = expected_data['locations']
#     #         compare_location_details(
#     #             locations, expected_locations,
#     #             ['start', 'end', 'representative'],
#     #             sig_acc, member_db
#     #         )

#     # for sig_acc, match_data in expected_dict.items():
#     #     assert sig_acc in result_dict, \
#     #         f"Signature {sig_acc} in the expected results but not actual results, {member_db}"

#     tree = ET.ElementTree(x_protein_elm[1])
#     ET.indent(tree, space="\t", level=0)
#     tree.write((x_matches_output_dir / f"{member_db}.matches.xml"), encoding='utf-8')


def test_cdd_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "CDD"
    match_data = load_match_data(member_db, x_matches_input_dir)

    result = xml_output.add_xml_output_matches(x_protein_elm[0], match_data)
    result_dict = xml_to_dict(result)
    result_dict = parse_matches_dict(result_dict['matches'][0])
    expected_dict = load_expected_match_data(member_db, x_matches_output_dir)

    for sig_acc, match_data in result_dict.items():
        assert sig_acc in expected_dict, f"Signature {sig_acc} not in expected results, {member_db}"
        if sig_acc in expected_dict:
            expected_data = expected_dict[sig_acc]
            compare_signature_details(match_data, expected_data, sig_acc, member_db)
            locations = match_data['locations']
            expected_locations = expected_data['locations']
            compare_location_details(
                locations, expected_locations,
                ['representative', 'evalue', 'score'],
                sig_acc, member_db
            )

    for sig_acc, match_data in expected_dict.items():
        assert sig_acc in result_dict, \
            f"Signature {sig_acc} in the expected results but not actual results, {member_db}"


def test_coils_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "COILS"
    match_data = load_match_data(member_db, x_matches_input_dir)

    result = xml_output.add_xml_output_matches(x_protein_elm[0], match_data)
    result_dict = xml_to_dict(result)
    result_dict = parse_matches_dict(result_dict['matches'][0])
    expected_dict = load_expected_match_data(member_db, x_matches_output_dir)

    for sig_acc, match_data in result_dict.items():
        assert sig_acc in expected_dict, f"Signature {sig_acc} not in expected results, {member_db}"
        if sig_acc in expected_dict:
            expected_data = expected_dict[sig_acc]
            compare_signature_details(match_data, expected_data, sig_acc, member_db)
            locations = match_data['locations']
            expected_locations = expected_data['locations']
            compare_location_details(
                locations, expected_locations,
                ['representative'],
                sig_acc, member_db
            )

    for sig_acc, match_data in expected_dict.items():
        assert sig_acc in result_dict, \
            f"Signature {sig_acc} in the expected results but not actual results, {member_db}"


def test_funfam_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "FUNFAM"
    match_data = load_match_data(member_db, x_matches_input_dir)

    result = xml_output.add_xml_output_matches(x_protein_elm[0], match_data)
    result_dict = xml_to_dict(result)
    result_dict = parse_matches_dict(result_dict['matches'][0])
    expected_dict = load_expected_match_data(member_db, x_matches_output_dir)

    for sig_acc, match_data in result_dict.items():
        assert sig_acc in expected_dict, f"Signature {sig_acc} not in expected results, {member_db}"
        if sig_acc in expected_dict:
            expected_data = expected_dict[sig_acc]
            compare_signature_details(match_data, expected_data, sig_acc, member_db)
            locations = match_data['locations']
            expected_locations = expected_data['locations']
            compare_location_details(
                locations, expected_locations,
                [
                    'representative', 'evalue', 'score',
                    'hmm-start', 'hmm-end', 'hmm-length', 'hmm-bounds',
                    'env-start', 'env-end', 'alignment', 'cigar-alignment'
                ],
                sig_acc, member_db
            )

    for sig_acc, match_data in expected_dict.items():
        assert sig_acc in result_dict, \
            f"Signature {sig_acc} in the expected results but not actual results, {member_db}"


def test_gene3d_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "GENE3D"
    match_data = load_match_data(member_db, x_matches_input_dir)

    result = xml_output.add_xml_output_matches(x_protein_elm[0], match_data)
    result_dict = xml_to_dict(result)
    result_dict = parse_matches_dict(result_dict['matches'][0])
    expected_dict = load_expected_match_data(member_db, x_matches_output_dir)

    for sig_acc, match_data in result_dict.items():
        assert sig_acc in expected_dict, f"Signature {sig_acc} not in expected results, {member_db}"
        if sig_acc in expected_dict:
            expected_data = expected_dict[sig_acc]
            compare_signature_details(match_data, expected_data, sig_acc, member_db)
            locations = match_data['locations']
            expected_locations = expected_data['locations']
            compare_location_details(
                locations, expected_locations,
                [
                    'representative', 'evalue', 'score',
                    'hmm-start', 'hmm-end', 'hmm-length', 'hmm-bounds',
                    'env-start', 'env-end', 'alignment', 'cigar-alignment'
                ],
                sig_acc, member_db
            )

    for sig_acc, match_data in expected_dict.items():
        assert sig_acc in result_dict, \
            f"Signature {sig_acc} in the expected results but not actual results, {member_db}"


def test_hamap_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "HAMAP"
    match_data = load_match_data(member_db, x_matches_input_dir)

    result = xml_output.add_xml_output_matches(x_protein_elm[0], match_data)
    result_dict = xml_to_dict(result)
    result_dict = parse_matches_dict(result_dict['matches'][0])
    expected_dict = load_expected_match_data(member_db, x_matches_output_dir)

    for sig_acc, match_data in result_dict.items():
        assert sig_acc in expected_dict, f"Signature {sig_acc} not in expected results, {member_db}"
        if sig_acc in expected_dict:
            expected_data = expected_dict[sig_acc]
            compare_signature_details(match_data, expected_data, sig_acc, member_db)
            locations = match_data['locations']
            expected_locations = expected_data['locations']
            assert locations[0]['score'] == expected_locations[0]['score'], \
                f"Mismatched score for {sig_acc} location, {member_db}"
            assert locations[0]['alignment'] == expected_locations[0]['alignment'], \
                f"Mismatched alignment for {sig_acc} location, {member_db}"

    for sig_acc, match_data in expected_dict.items():
        assert sig_acc in result_dict, \
            f"Signature {sig_acc} in the expected results but not actual results, {member_db}"


def test_mobidb_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "MOBIDB"
    match_data = load_match_data(member_db, x_matches_input_dir)

    result = xml_output.add_xml_output_matches(x_protein_elm[0], match_data)
    result_dict = xml_to_dict(result)
    result_dict = parse_matches_dict(result_dict['matches'][0])
    expected_dict = load_expected_match_data(member_db, x_matches_output_dir)

    for sig_acc, match_data in result_dict.items():
        assert sig_acc in expected_dict, f"Signature {sig_acc} not in expected results, {member_db}"
        if sig_acc in expected_dict:
            expected_data = expected_dict[sig_acc]
            compare_signature_details(match_data, expected_data, sig_acc, member_db)
            locations = match_data['locations']
            expected_locations = expected_data['locations']
            compare_location_details(
                locations, expected_locations,
                ['start', 'end', 'representative'],
                sig_acc, member_db
            )

    for sig_acc, match_data in expected_dict.items():
        assert sig_acc in result_dict, \
            f"Signature {sig_acc} in the expected results but not actual results, {member_db}"


def test_ncbifam_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "NCBIFAM"
    match_data = load_match_data(member_db, x_matches_input_dir)

    result = xml_output.add_xml_output_matches(x_protein_elm[0], match_data)
    result_dict = xml_to_dict(result)
    result_dict = parse_matches_dict(result_dict['matches'][0])
    expected_dict = load_expected_match_data(member_db, x_matches_output_dir)

    for sig_acc, match_data in result_dict.items():
        assert sig_acc in expected_dict, f"Signature {sig_acc} not in expected results, {member_db}"
        if sig_acc in expected_dict:
            expected_data = expected_dict[sig_acc]
            compare_signature_details(match_data, expected_data, sig_acc, member_db)
            locations = match_data['locations']
            expected_locations = expected_data['locations']
            compare_location_details(
                locations, expected_locations,
                [
                    'start', 'end', 'representative',
                    'env-start', 'env-end', 'score', 'evalue',
                    'hmm-start', 'hmm-end', 'hmm-length',
                    'hmm-bounds', 'alignment', 'cigar-alignment'
                ],
                sig_acc, member_db
            )

    for sig_acc, match_data in expected_dict.items():
        assert sig_acc in result_dict, \
            f"Signature {sig_acc} in the expected results but not actual results, {member_db}"


def test_panther_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "PANTHER"
    match_data = load_match_data(member_db, x_matches_input_dir)

    result = xml_output.add_xml_output_matches(x_protein_elm[0], match_data)
    result_dict = xml_to_dict(result)
    result_dict = parse_matches_dict(result_dict['matches'][0])
    expected_dict = load_expected_match_data(member_db, x_matches_output_dir)

    for sig_acc, match_data in result_dict.items():
        assert sig_acc in expected_dict, f"Signature {sig_acc} not in expected results, {member_db}"
        if sig_acc in expected_dict:
            expected_data = expected_dict[sig_acc]
            compare_signature_details(match_data, expected_data, sig_acc, member_db)
            locations = match_data['locations']
            expected_locations = expected_data['locations']
            compare_location_details(
                locations, expected_locations,
                [
                    'start', 'end', 'representative',
                    'env-start', 'env-end', 'score', 'evalue',
                    'hmm-start', 'hmm-end', 'hmm-length',
                    'hmm-bounds', 'alignment', 'cigar-alignment'
                ],
                sig_acc, member_db
            )

    for sig_acc, match_data in expected_dict.items():
        assert sig_acc in result_dict, \
            f"Signature {sig_acc} in the expected results but not actual results, {member_db}"

    tree = ET.ElementTree(x_protein_elm[1])
    ET.indent(tree, space="\t", level=0)
    tree.write((x_matches_output_dir / f"{member_db}.matches.xml"), encoding='utf-8')


def test_pfam_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "PFAM"
    match_data = load_match_data(member_db, x_matches_input_dir)

    result = xml_output.add_xml_output_matches(x_protein_elm[0], match_data)
    result_dict = xml_to_dict(result)
    result_dict = parse_matches_dict(result_dict['matches'][0])
    expected_dict = load_expected_match_data(member_db, x_matches_output_dir)

    for sig_acc, match_data in result_dict.items():
        assert sig_acc in expected_dict, f"Signature {sig_acc} not in expected results, {member_db}"
        if sig_acc in expected_dict:
            expected_data = expected_dict[sig_acc]
            compare_signature_details(match_data, expected_data, sig_acc, member_db)
            locations = match_data['locations']
            expected_locations = expected_data['locations']
            compare_location_details(
                locations, expected_locations,
                ['start', 'end', 'representative'],
                sig_acc, member_db
            )

    for sig_acc, match_data in expected_dict.items():
        assert sig_acc in result_dict, \
            f"Signature {sig_acc} in the expected results but not actual results, {member_db}"


def test_phobius_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "PHOBIUS"
    match_data = load_match_data(member_db, x_matches_input_dir)

    result = xml_output.add_xml_output_matches(x_protein_elm[0], match_data)
    result_dict = xml_to_dict(result)
    result_dict = parse_matches_dict(result_dict['matches'][0])
    expected_dict = load_expected_match_data(member_db, x_matches_output_dir)

    for sig_acc, match_data in result_dict.items():
        assert sig_acc in expected_dict, f"Signature {sig_acc} not in expected results, {member_db}"
        if sig_acc in expected_dict:
            expected_data = expected_dict[sig_acc]
            compare_signature_details(match_data, expected_data, sig_acc, member_db)
            locations = match_data['locations']
            expected_locations = expected_data['locations']
            compare_location_details(
                locations, expected_locations,
                ['start', 'end', 'representative'],
                sig_acc, member_db
            )

    for sig_acc, match_data in expected_dict.items():
        assert sig_acc in result_dict, \
            f"Signature {sig_acc} in the expected results but not actual results, {member_db}"


def test_pirsf_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "PIRSF"
    match_data = load_match_data(member_db, x_matches_input_dir)

    result = xml_output.add_xml_output_matches(x_protein_elm[0], match_data)
    result_dict = xml_to_dict(result)
    result_dict = parse_matches_dict(result_dict['matches'][0])
    expected_dict = load_expected_match_data(member_db, x_matches_output_dir)

    for sig_acc, match_data in result_dict.items():
        assert sig_acc in expected_dict, f"Signature {sig_acc} not in expected results, {member_db}"
        if sig_acc in expected_dict:
            expected_data = expected_dict[sig_acc]
            compare_signature_details(match_data, expected_data, sig_acc, member_db)
            locations = match_data['locations']
            expected_locations = expected_data['locations']
            compare_location_details(
                locations, expected_locations,
                [
                    'start', 'end', 'representative',
                    'env-end', 'env-start', 'evalue', 'score',
                    'hmm-start', 'hmm-end', 'hmm-length', 'hmm-bounds'
                ],
                sig_acc, member_db
            )

    for sig_acc, match_data in expected_dict.items():
        assert sig_acc in result_dict, \
            f"Signature {sig_acc} in the expected results but not actual results, {member_db}"


# def test_pirsr_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
#     member_db = "PIRSR"
#     match_data = load_match_data(member_db, x_matches_input_dir)

#     result = xml_output.add_xml_output_matches(x_protein_elm[0], match_data)
#     result_dict = xml_to_dict(result)
#     result_dict = parse_matches_dict(result_dict['matches'][0])
#     expected_dict = load_expected_match_data(member_db, x_matches_output_dir)

#     for sig_acc, match_data in result_dict.items():
#         assert sig_acc in expected_dict, f"Signature {sig_acc} not in expected results, {member_db}"
#         if sig_acc in expected_dict:
#             expected_data = expected_dict[sig_acc]
#             compare_signature_details(match_data, expected_data, sig_acc, member_db)
#             locations = match_data['locations']
#             expected_locations = expected_data['locations']
#             compare_location_details(
#                 locations, expected_locations,
#                 [
#                     'start', 'end', 'representative',
#                     'env-end', 'env-start', 'evalue', 'score',
#                     'hmm-start', 'hmm-end', 'hmm-length', 'hmm-bounds'
#                 ],
#                 sig_acc, member_db
#             )

#     for sig_acc, match_data in expected_dict.items():
#         assert sig_acc in result_dict, \
#             f"Signature {sig_acc} in the expected results but not actual results, {member_db}"

#     tree = ET.ElementTree(x_protein_elm[1])
#     ET.indent(tree, space="\t", level=0)
#     tree.write((x_matches_output_dir / f"{member_db}.matches.xml"), encoding='utf-8')


def test_prints_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "PRINTS"
    match_data = load_match_data(member_db, x_matches_input_dir)

    result = xml_output.add_xml_output_matches(x_protein_elm[0], match_data)
    result_dict = xml_to_dict(result)
    result_dict = parse_matches_dict(result_dict['matches'][0])
    expected_dict = load_expected_match_data(member_db, x_matches_output_dir)

    for sig_acc, match_data in result_dict.items():
        assert sig_acc in expected_dict, f"Signature {sig_acc} not in expected results, {member_db}"
        if sig_acc in expected_dict:
            expected_data = expected_dict[sig_acc]
            compare_signature_details(match_data, expected_data, sig_acc, member_db)
            locations = match_data['locations']
            expected_locations = expected_data['locations']
            compare_location_details(
                locations, expected_locations,
                [
                    'motifNumber', 'pvalue', 'score',
                    'start', 'end', 'representative'
                ],
                sig_acc, member_db
            )

    for sig_acc, match_data in expected_dict.items():
        assert sig_acc in result_dict, \
            f"Signature {sig_acc} in the expected results but not actual results, {member_db}"


def test_prosite_pattern_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "PROSITE_PATTERNS"
    match_data = load_match_data(member_db, x_matches_input_dir)

    result = xml_output.add_xml_output_matches(x_protein_elm[0], match_data)
    result_dict = xml_to_dict(result)
    result_dict = parse_matches_dict(result_dict['matches'][0])
    expected_dict = load_expected_match_data(member_db, x_matches_output_dir)

    for sig_acc, match_data in result_dict.items():
        assert sig_acc in expected_dict, f"Signature {sig_acc} not in expected results, {member_db}"
        if sig_acc in expected_dict:
            expected_data = expected_dict[sig_acc]
            compare_signature_details(match_data, expected_data, sig_acc, member_db)
            locations = match_data['locations']
            expected_locations = expected_data['locations']
            compare_location_details(
                locations, expected_locations,
                ['start', 'end', 'alignment', 'cigar-alignment'],
                sig_acc, member_db
            )

    for sig_acc, match_data in expected_dict.items():
        assert sig_acc in result_dict, \
            f"Signature {sig_acc} in the expected results but not actual results, {member_db}"


def test_prosite_profile_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "PROSITE_PROFILES"
    match_data = load_match_data(member_db, x_matches_input_dir)

    result = xml_output.add_xml_output_matches(x_protein_elm[0], match_data)
    result_dict = xml_to_dict(result)
    result_dict = parse_matches_dict(result_dict['matches'][0])
    expected_dict = load_expected_match_data(member_db, x_matches_output_dir)

    for sig_acc, match_data in result_dict.items():
        assert sig_acc in expected_dict, f"Signature {sig_acc} not in expected results, {member_db}"
        if sig_acc in expected_dict:
            expected_data = expected_dict[sig_acc]
            compare_signature_details(match_data, expected_data, sig_acc, member_db)
            locations = match_data['locations']
            expected_locations = expected_data['locations']
            compare_location_details(
                locations, expected_locations,
                ['alignment', 'score', 'start', 'end'],
                sig_acc, member_db
            )

    for sig_acc, match_data in expected_dict.items():
        assert sig_acc in result_dict, \
            f"Signature {sig_acc} in the expected results but not actual results, {member_db}"


# def test_sfld_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
#     member_db = "SFLD"
#     match_data = load_match_data(member_db, x_matches_input_dir)

#     result = xml_output.add_xml_output_matches(x_protein_elm[0], match_data)
#     result_dict = xml_to_dict(result)
#     result_dict = parse_matches_dict(result_dict['matches'][0])
#     expected_dict = load_expected_match_data(member_db, x_matches_output_dir)

#     for sig_acc, match_data in result_dict.items():
#         assert sig_acc in expected_dict, f"Signature {sig_acc} not in expected results, {member_db}"
#         if sig_acc in expected_dict:
#             expected_data = expected_dict[sig_acc]
#             compare_signature_details(match_data, expected_data, sig_acc, member_db)
#             locations = match_data['locations']
#             expected_locations = expected_data['locations']
#             compare_location_details(
#                 locations, expected_locations,
#                 ['alignment', 'score', 'start', 'end'],
#                 sig_acc, member_db
#             )

#     for sig_acc, match_data in expected_dict.items():
#         assert sig_acc in result_dict, \
#             f"Signature {sig_acc} in the expected results but not actual results, {member_db}"


def test_signalp_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "SIGNALP"
    match_data = load_match_data(member_db, x_matches_input_dir)

    result = xml_output.add_xml_output_matches(x_protein_elm[0], match_data)
    result_dict = xml_to_dict(result)
    result_dict = parse_matches_dict(result_dict['matches'][0])
    expected_dict = load_expected_match_data(member_db, x_matches_output_dir)

    for sig_acc, match_data in result_dict.items():
        assert sig_acc in expected_dict, f"Signature {sig_acc} not in expected results, {member_db}"
        if sig_acc in expected_dict:
            expected_data = expected_dict[sig_acc]
            compare_signature_details(match_data, expected_data, sig_acc, member_db)
            locations = match_data['locations']
            expected_locations = expected_data['locations']
            compare_location_details(
                locations, expected_locations,
                [
                    'end', 'start', 'pvalue', 'cleavage_start',
                    'cleavage_end', 'representative'
                ],
                sig_acc, member_db
            )

    for sig_acc, match_data in expected_dict.items():
        assert sig_acc in result_dict, \
            f"Signature {sig_acc} in the expected results but not actual results, {member_db}"


def test_smart_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "SMART"
    match_data = load_match_data(member_db, x_matches_input_dir)

    result = xml_output.add_xml_output_matches(x_protein_elm[0], match_data)
    result_dict = xml_to_dict(result)
    result_dict = parse_matches_dict(result_dict['matches'][0])
    expected_dict = load_expected_match_data(member_db, x_matches_output_dir)

    for sig_acc, match_data in result_dict.items():
        assert sig_acc in expected_dict, f"Signature {sig_acc} not in expected results, {member_db}"
        if sig_acc in expected_dict:
            expected_data = expected_dict[sig_acc]
            compare_signature_details(match_data, expected_data, sig_acc, member_db)
            locations = match_data['locations']
            expected_locations = expected_data['locations']
            compare_location_details(
                locations, expected_locations,
                [
                    'score', 'evalue', 'hmm-start', 'hmm-end',
                    'hmm-len', 'hmm-bounds',
                    'start', 'end', 'representative'
                ],
                sig_acc, member_db
            )

    for sig_acc, match_data in expected_dict.items():
        assert sig_acc in result_dict, \
            f"Signature {sig_acc} in the expected results but not actual results, {member_db}"


def test_superfamily_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "SUPERFAMILY"
    match_data = load_match_data(member_db, x_matches_input_dir)

    result = xml_output.add_xml_output_matches(x_protein_elm[0], match_data)
    result_dict = xml_to_dict(result)
    result_dict = parse_matches_dict(result_dict['matches'][0])
    expected_dict = load_expected_match_data(member_db, x_matches_output_dir)

    for sig_acc, match_data in result_dict.items():
        assert sig_acc in expected_dict, f"Signature {sig_acc} not in expected results, {member_db}"
        if sig_acc in expected_dict:
            expected_data = expected_dict[sig_acc]
            compare_signature_details(match_data, expected_data, sig_acc, member_db)
            locations = match_data['locations']
            expected_locations = expected_data['locations']
            compare_location_details(
                locations, expected_locations,
                ['start', 'end', 'representative'],
                sig_acc, member_db
            )

    for sig_acc, match_data in expected_dict.items():
        assert sig_acc in result_dict, \
            f"Signature {sig_acc} in the expected results but not actual results, {member_db}"
