"""Test the python script that writes the XML output.

Fixtures are prefixed with 'x_' to indicate they are used in unit
testing the XML format writer.

These test are intened to be run from the root of the repository using:
python -m pytest -v
"""

import json
import re

import xml.etree.ElementTree as ET

import pytest

from interproscan.scripts.output.format_writer import xml_output


PIRSR_SITE_LOCATION = re.compile(r"'end':\s(\d+),\s'residue':\s'(\D|\[\D+\])',\s'start':\s(\d+)")
PIRSR_SITES_REGEX = re.compile(r"'description':\s'(.*?)',\s'group':\s(\d+),\s'hmmEnd':\s(\d+),\s'hmmStart':\s(\d+),\s'label':\s(.*?),\s'numLocations':\s(\d+),\s'siteLocations':\s(\[\{.*?\}\])")
SFLD_SITE_LOCATION = re.compile(r"\{'start':\s(\d+),\s'end':\s(\d+),\s'residue':\s'(\D)'\}")
SFLD_SITES_REGEX = re.compile(r"'description':\s'(.*?)',\s'group':\s(\d+),\s'hmmEnd':\s(\d+),\s'hmmStart':\s(\d+),\s'label':\s(.*?),\s'numLocations':\s(\d+),\s'siteLocations':\s(\[\{.*?'\}\])")


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
                if 'analysis-location' in location:
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


def compare_xml_matches(
    member_db: str,
    matches_input_dir: pytest.fixture,
    matches_output_dir: pytest.fixture,
    protein_elm: pytest.fixture,
    location_keys: list[str],
) -> None:
    match_data = load_match_data(member_db, matches_input_dir)
    result = xml_output.add_xml_output_matches(protein_elm[0], match_data)
    result_dict = xml_to_dict(result)
    result_dict = parse_matches_dict(result_dict['matches'][0])
    expected_dict = load_expected_match_data(member_db, matches_output_dir)
    for sig_acc, match_data in result_dict.items():
        assert sig_acc in expected_dict, f"Signature {sig_acc} not in expected results, {member_db}"
        if sig_acc in expected_dict:
            expected_data = expected_dict[sig_acc]
            compare_signature_details(match_data, expected_data, sig_acc, member_db)
            locations = match_data['locations']
            expected_locations = expected_data['locations']
            compare_location_details(
                locations, expected_locations,
                location_keys,
                sig_acc, member_db
            )

    for sig_acc, match_data in expected_dict.items():
        assert sig_acc in result_dict, \
            f"Signature {sig_acc} in the expected results but not actual results, {member_db}"


def test_antifam_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "ANTIFAM"
    compare_xml_matches(
        member_db,
        x_matches_input_dir,
        x_matches_output_dir,
        x_protein_elm,
        ['start', 'end', 'representative']
    )


def test_cdd_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "CDD"
    compare_xml_matches(
        member_db,
        x_matches_input_dir,
        x_matches_output_dir,
        x_protein_elm,
        ['representative', 'evalue', 'score']
    )


def test_coils_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "COILS"
    compare_xml_matches(
        member_db,
        x_matches_input_dir,
        x_matches_output_dir,
        x_protein_elm,
        ['representative']
    )


def test_funfam_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "FUNFAM"
    compare_xml_matches(
        member_db,
        x_matches_input_dir,
        x_matches_output_dir,
        x_protein_elm,
        [
            'representative', 'evalue', 'score',
            'hmm-start', 'hmm-end', 'hmm-length', 'hmm-bounds',
            'env-start', 'env-end', 'alignment', 'cigar-alignment'
        ]
    )


def test_gene3d_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "GENE3D"
    compare_xml_matches(
        member_db,
        x_matches_input_dir,
        x_matches_output_dir,
        x_protein_elm,
        [
            'representative', 'evalue', 'score',
            'hmm-start', 'hmm-end', 'hmm-length', 'hmm-bounds',
            'env-start', 'env-end', 'alignment', 'cigar-alignment'
        ]
    )


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
    member_db = "MOBIDB_LITE"
    compare_xml_matches(
        member_db,
        x_matches_input_dir,
        x_matches_output_dir,
        x_protein_elm,
        ['start', 'end', 'representative']
    )


def test_ncbifam_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "NCBIFAM"
    compare_xml_matches(
        member_db,
        x_matches_input_dir,
        x_matches_output_dir,
        x_protein_elm,
        [
            'start', 'end', 'representative',
            'env-start', 'env-end', 'score', 'evalue',
            'hmm-start', 'hmm-end', 'hmm-length',
            'hmm-bounds', 'alignment', 'cigar-alignment'
        ]
    )


def test_panther_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "PANTHER"
    compare_xml_matches(
        member_db,
        x_matches_input_dir,
        x_matches_output_dir,
        x_protein_elm,
        [
            'start', 'end', 'representative',
            'env-start', 'env-end', 'score', 'evalue',
            'hmm-start', 'hmm-end', 'hmm-length',
            'hmm-bounds', 'alignment', 'cigar-alignment'
        ]
    )


def test_pfam_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "PFAM"
    compare_xml_matches(
        member_db,
        x_matches_input_dir,
        x_matches_output_dir,
        x_protein_elm,
        ['start', 'end', 'representative']
    )


def test_phobius_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "PHOBIUS"
    compare_xml_matches(
        member_db,
        x_matches_input_dir,
        x_matches_output_dir,
        x_protein_elm,
        ['start', 'end', 'representative']
    )


def test_pirsf_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "PIRSF"
    compare_xml_matches(
        member_db,
        x_matches_input_dir,
        x_matches_output_dir,
        x_protein_elm,
        [
            'start', 'end', 'representative',
            'env-end', 'env-start', 'evalue', 'score',
            'hmm-start', 'hmm-end', 'hmm-length', 'hmm-bounds'
        ]
    )


def test_pirsr_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "PIRSR"
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

            for location in locations:
                if 'sites' not in location:
                    continue
                for rsult_site in location['sites']:
                    site_found = False
                    # rsult_site is a dict {'description':.., 'group':.., 'hmmEnd':....}
                    
                    for xpctd_location in expected_locations:
                        # xpcted_location is a str rerpr of the rsult_site dict
                        expected_sites = PIRSR_SITES_REGEX.findall(xpctd_location['sites'])
                        for expected_site in expected_sites:
                            # expected_site is a tuple
                            if int(rsult_site['hmmStart']) == int(expected_site[3]) and \
                                int(rsult_site['hmmEnd']) == int(expected_site[2]) and \
                                rsult_site['description'].strip() == expected_site[0].strip():
                                site_found = True
                                assert str(rsult_site['group']) == expected_site[1], \
                                    (
                                        "Mismatch site 'group' value for site hmmStart:"
                                        f"{rsult_site['hmmStart']}, hmmEnd:{rsult_site['hmmEnd']},"
                                        f"{sig_acc}, {member_db}"
                                    )
                                assert str(rsult_site['label']) == expected_site[4].strip("'"), \
                                    (
                                        "Mismatch site 'label' value for site hmmStart:"
                                        f"{rsult_site['hmmStart']}, hmmEnd:{rsult_site['hmmEnd']},"
                                        f"{sig_acc}, {member_db}"
                                    )
                                assert rsult_site['numLocations'] == int(expected_site[5]), \
                                    (
                                        "Mismatch site number of locations for site hmmStart:"
                                        f"{rsult_site['hmmStart']}, hmmEnd:{rsult_site['hmmEnd']},"
                                        f"{sig_acc}, {member_db}"
                                    )

                                # check the site locations
                                expected_site_locations = PIRSR_SITE_LOCATION.findall(expected_site[-1])
                                for r_sloc in rsult_site['siteLocations']:
                                    site_loc_found = False
                                    for ex_sloc in expected_site_locations:
                                        if str(r_sloc['start']) == str(ex_sloc[2]) and \
                                            str(r_sloc['end']) == str(ex_sloc[0]) and \
                                            r_sloc['residue'] == ex_sloc[1]:
                                            site_loc_found = True
                                    if not site_loc_found:
                                        pytest.fail(
                                            "Could not find site location in current results ("
                                            f"'{r_sloc['start']}', '{r_sloc['end']}', "
                                            f"'{r_sloc['residue']}') in expected results"
                                        )
                                for ex_sloc in expected_site_locations:
                                    site_loc_found = False
                                    for r_sloc in rsult_site['siteLocations']:
                                        if str(r_sloc['start']) == str(ex_sloc[2]) and \
                                            str(r_sloc['end']) == str(ex_sloc[0]) and \
                                            r_sloc['residue'] == ex_sloc[1]:
                                            site_loc_found = True
                                    if not site_loc_found:
                                        pytest.fail(
                                            "Could not find site location in current results ("
                                            f"'{r_sloc['start']}', '{r_sloc['end']}', "
                                            f"'{r_sloc['residue']}') in expected results"
                                        )
                    if not site_found:
                        pytest.fail(
                            f"Result site '{rsult_site['description']}' with "
                            f"hmmStart {rsult_site['hmmStart']} and hmmEnd "
                            f"{rsult_site['hmmEnd']} for {sig_acc}, {member_db}, "
                            "not found in expected results"
                        )

    for sig_acc, match_data in expected_dict.items():
        assert sig_acc in result_dict, \
            f"Signature {sig_acc} in the expected results but not actual results, {member_db}"


def test_prints_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "PRINTS"
    compare_xml_matches(
        member_db,
        x_matches_input_dir,
        x_matches_output_dir,
        x_protein_elm,
        [
            'motifNumber', 'pvalue', 'score',
            'start', 'end', 'representative'
        ]
    )


def test_prosite_pattern_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "PROSITE_PATTERNS"
    compare_xml_matches(
        member_db,
        x_matches_input_dir,
        x_matches_output_dir,
        x_protein_elm,
        ['start', 'end', 'alignment', 'cigar-alignment']
    )


def test_prosite_profile_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "PROSITE_PROFILES"
    compare_xml_matches(
        member_db,
        x_matches_input_dir,
        x_matches_output_dir,
        x_protein_elm,
        ['alignment', 'score', 'start', 'end']
    )


def test_sfld_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "SFLD"
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

            for location in locations:
                for rsult_site in location['sites']:
                    site_found = False
                    # rsult_site is a dict {'description':.., 'group':.., 'hmmEnd':....}
                    
                    for xpctd_location in expected_locations:
                        # xpcted_location is a str rerpr of the rsult_site dict
                        expected_sites = SFLD_SITES_REGEX.findall(xpctd_location['sites'])
                        for expected_site in expected_sites:
                            # expected_site is a tuple
                            if int(rsult_site['hmmStart']) == int(expected_site[3]) and \
                                int(rsult_site['hmmEnd']) == int(expected_site[2]) and \
                                rsult_site['description'] == expected_site[0]:
                                site_found = True
                                assert str(rsult_site['group']) == expected_site[1], \
                                    (
                                        "Mismatch site 'group' value for site hmmStart:"
                                        f"{rsult_site['hmmStart']}, hmmEnd:{rsult_site['hmmEnd']},"
                                        f"{sig_acc}, {member_db}"
                                    )
                                assert str(rsult_site['label']) == expected_site[4], \
                                    (
                                        "Mismatch site 'label' value for site hmmStart:"
                                        f"{rsult_site['hmmStart']}, hmmEnd:{rsult_site['hmmEnd']},"
                                        f"{sig_acc}, {member_db}"
                                    )
                                assert rsult_site['numLocations'] == int(expected_site[5]), \
                                    (
                                        "Mismatch site number of locations for site hmmStart:"
                                        f"{rsult_site['hmmStart']}, hmmEnd:{rsult_site['hmmEnd']},"
                                        f"{sig_acc}, {member_db}"
                                    )
                                
                                # check the site locations
                                expected_site_locations = SFLD_SITE_LOCATION.findall(expected_site[-1])
                                for r_sloc in rsult_site['siteLocations']:
                                    site_loc_found = False
                                    for ex_sloc in expected_site_locations:
                                        if r_sloc['start'] == int(ex_sloc[0]) and \
                                            r_sloc['end'] == int(ex_sloc[1]) and \
                                            r_sloc['residue'] == ex_sloc[2]:
                                            site_loc_found = True
                                    if not site_loc_found:
                                        pytest.fail(
                                            "Could not find site location in current results ("
                                            f"'{r_sloc['start']}', '{r_sloc['end']}', "
                                            f"'{r_sloc['residue']}') in expected results"
                                        )
                                for ex_sloc in expected_site_locations:
                                    site_loc_found = False
                                    for r_sloc in rsult_site['siteLocations']:
                                        if r_sloc['start'] == int(ex_sloc[0]) and \
                                            r_sloc['end'] == int(ex_sloc[1]) and \
                                            r_sloc['residue'] == ex_sloc[2]:
                                            site_loc_found = True
                                    if not site_loc_found:
                                        pytest.fail(
                                            "Could not find site location in current results ("
                                            f"'{r_sloc['start']}', '{r_sloc['end']}', "
                                            f"'{r_sloc['residue']}') in expected results"
                                        )
                    if not site_found:
                        pytest.fail(
                            f"Result site with hmmStart {rsult_site['hmmStart']} and hmmEnd "
                            f"{rsult_site['hmmEnd']} for {sig_acc}, {member_db}, "
                            "not found in expected results"
                        )

    for sig_acc, match_data in expected_dict.items():
        assert sig_acc in result_dict, \
            f"Signature {sig_acc} in the expected results but not actual results, {member_db}"


def test_signalp_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "SIGNALP"
    compare_xml_matches(
        member_db,
        x_matches_input_dir,
        x_matches_output_dir,
        x_protein_elm,
        [
            'end', 'start', 'pvalue', 'cleavage_start',
            'cleavage_end', 'representative'
        ]
    )


def test_smart_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "SMART"
    compare_xml_matches(
        member_db,
        x_matches_input_dir,
        x_matches_output_dir,
        x_protein_elm,
        [
            'score', 'evalue', 'hmm-start', 'hmm-end',
            'hmm-len', 'hmm-bounds',
            'start', 'end', 'representative'
        ]
    )


def test_superfamily_xml_match(x_matches_input_dir, x_matches_output_dir, x_protein_elm):
    member_db = "SUPERFAMILY"
    compare_xml_matches(
        member_db,
        x_matches_input_dir,
        x_matches_output_dir,
        x_protein_elm,
        ['start', 'end', 'representative']
    )
