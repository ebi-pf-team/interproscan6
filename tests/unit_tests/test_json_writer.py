"""Test the python script that coordinates writing the output.

Fixtures are prefixed with 'j_' to indicate they are used in unit
testing the JSON format writer.

These test are intened to be run from the root of the repository using:
python -m pytest -v
"""

import json

import pytest

from interproscan.scripts.output.format_writer import json_output


@pytest.fixture
def j_expected_nucleic(test_output_dir):
    _path = test_output_dir / "format_writer/build_funcs/expected_json_nucleic.json"
    with open(_path, "r") as fh:
        expected_output = json.load(fh)
    return expected_output

@pytest.fixture
def j_expected_protein(test_output_dir):
    _path = test_output_dir / "format_writer/build_funcs/expected_json_protein.json"
    with open(_path, "r") as fh:
        expected_output = json.load(fh)
    return expected_output


@pytest.fixture
def j_matches_input_dir(test_input_dir):
    return test_input_dir / "format_writer/matches"


@pytest.fixture
def j_matches_output_dir(test_output_dir):
    return test_output_dir / "format_writer/matches/json/"


@pytest.fixture
def j_nucleic_seq_matches(test_input_dir):
    _path = test_input_dir / "format_writer/nucleic/seq_match_dict.json"
    with open(_path, "r") as fh:
        _input = json.load(fh)
    return _input


@pytest.fixture
def j_out_path(test_output_dir):
    outdir = test_output_dir / "temp"
    outdir.mkdir(parents=True, exist_ok=True)
    filepath = outdir / "tsv.unittest.ips6.tsv"
    return filepath


@pytest.fixture
def j_prot_seq_matches(test_input_dir):
    _path = test_input_dir / "format_writer/protein/seq_match_dict.json"
    with open(_path, "r") as fh:
        _input = json.load(fh)
    return _input


def load_expected_match_data(member_db, matches_output_dir):
    with open((matches_output_dir / f"{member_db.upper()}.match.json"), "r") as fh:
        match_list = json.load(fh)
    return match_list


def load_match_data(member_db, matches_input_dir):
    with open((matches_input_dir / f"{member_db.upper()}.match-data.json"), "r") as fh:
        match_data = json.load(fh)
    return match_data


def parse_matches_into_dict(match_list: list[dict[str, dict]]) -> dict[str, dict[str, dict]]:
    match_dict = {}
    for match in match_list:
        sig_acc = match["signature"]["accession"]
        match_dict[sig_acc] = match
    return match_dict


def test_build_json_protein(j_prot_seq_matches, j_out_path, j_expected_protein, monkeypatch):
    """Tests the function runs.
    The output of get_matches() is assessed in its own function and the overall
    output is assessed in the integration tests.
    """
    def mock_get_matches(*args, **kwards):
        return []

    monkeypatch.setattr(json_output, "get_matches", mock_get_matches)

    json_output.build_json_output_protein(j_prot_seq_matches, j_out_path, "unit.test.version")

    with open(j_out_path, "r") as fh:
        result = json.load(fh)

    assert result['interproscan-version'] == j_expected_protein['interproscan-version'], "InterProScan version does not match"
    assert len(result['results']) == len(j_expected_protein['results']), "Number of items in 'results' between the current and expected output do not match"
    # there aren't any matches to check for matching at this point because of the mocking

    j_out_path.unlink()


def test_build_json_nucleic(j_nucleic_seq_matches, j_out_path, j_expected_nucleic, monkeypatch):
    """Tests the function runs.
    The output of get_matches() is assessed in its own function and the overall
    output is assessed in the integration tests.
    """
    def mock_get_matches(*args, **kwards):
        return []

    monkeypatch.setattr(json_output, "get_matches", mock_get_matches)

    json_output.build_json_output_nucleic(j_nucleic_seq_matches, j_out_path, "unit.test.version")

    with open(j_out_path, "r") as fh:
        result = json.load(fh)

    assert result['interproscan-version'] == j_expected_nucleic['interproscan-version'], "InterProScan version does not match"
    assert len(result['results']) == len(j_expected_nucleic['results']), "Number of items in 'results' between the current and expected output do not match"
    # there aren't any matches to check for matching at this point because of the mocking

    j_out_path.unlink()


def compare_signature_details(match_data: dict, expected_data: dict, sig_acc: str, member_db: str) -> None:
    assert match_data['model-ac'] == expected_data['model-ac'], f"Mismatched 'model-ac for {sig_acc}, {member_db}"
    assert all(
        match_data['signature'][key] == expected_data['signature'][key]
        for key in ['accession', 'name', 'description']
    ), f"Mismatch in 'signature' details for {sig_acc}, {member_db}"
    assert all(
        match_data['signature']['signatureLibraryRelease'][key] == expected_data['signature']['signatureLibraryRelease'][key]
        for key in ['library', 'version']
    ), f"Mismatch in 'signatureLibraryRelease' details for {sig_acc}, {member_db}"



def test_cdd_json_match(j_matches_input_dir, j_matches_output_dir):
    member_db = "CDD"
    match_data = load_match_data(member_db, j_matches_input_dir)

    result_list = json_output.get_matches(match_data)

    expected_list = load_expected_match_data(member_db, j_matches_output_dir)

    # check the number of matches is the same
    assert len(result_list) == len(expected_list), "Different number of matches retrieved"

    # check the match data matches
    result_dict = parse_matches_into_dict(result_list)
    expected_dict = parse_matches_into_dict(expected_list)

    for sig_acc, match_data in result_dict. items():
        assert sig_acc in expected_dict, f"Signature {sig_acc} not in expected results"
        if sig_acc in expected_dict:
            # compare the signature details
            expected_data = expected_dict[sig_acc]
            compare_signature_details(match_data, expected_data, sig_acc, member_db)

            # compare locations
            assert len(match_data['locations']) == len(expected_data['locations']), f"Mismatched number of locations for {sig_acc}, {member_db}"
            for location in match_data['locations']:
                found = False
                for expected_location in expected_data['locations']:
                    if location['start'] == expected_location['start'] \
                        and location['end'] == expected_location['end']:
                        found = True
                        assert all(
                            location[key] == expected_location[key]
                            for key in ['representative', 'evalue', 'score']
                        ), (
                            f"Mismatch in 'location' details for {sig_acc}, {member_db} "
                            f"(start: {location['start']} - end {location['end']})"
                        )
                        assert len(location['sites']) == len(expected_location['sites'])
                        assert len(location['location-fragments']) == len(expected_location['location-fragments'])
                        # gritty detailed comparison can be handeled by the intergration test
                if not found:
                    pytest.fail((
                        f"Location with start {location['start']} and end {location['end']} "
                        f"not found in expected data for {sig_acc}, {member_db}"
                    ))
            for expected_location in expected_data['locations']:
                found = False
                for location in match_data['locations']:
                    if location['start'] == expected_location['start'] \
                        and location['end'] == expected_location['end']:
                        found = True
                if not found:
                    pytest.fail((
                        f"Expected location with start {location['start']} and end {location['end']} "
                        f"not found in test result data for {sig_acc}, {member_db}"
                    ))


# def test_coils_json_match(j_matches_input_dir):
#     match_data = load_match_data("COILS", j_matches_input_dir)

#     result = json_output.get_matches(match_data)

#     with open("tests/unit_tests/test_outputs/format_writer/matches/json/COILS.match.json", "w") as fh:
#         json.dump(result, fh)


# def test_funfam_json_match(j_matches_input_dir):
#     match_data = load_match_data("FUNFAM", j_matches_input_dir)

#     result = json_output.get_matches(match_data)

#     with open("tests/unit_tests/test_outputs/format_writer/matches/json/FUNFAM.match.json", "w") as fh:
#         json.dump(result, fh)


# def test_gene3d_json_match(j_matches_input_dir):
#     match_data = load_match_data("GENE3D", j_matches_input_dir)

#     result = json_output.get_matches(match_data)

#     with open("tests/unit_tests/test_outputs/format_writer/matches/json/GENE3D.match.json", "w") as fh:
#         json.dump(result, fh)


# def test_hamap_json_match(j_matches_input_dir):
#     match_data = load_match_data("HAMAP", j_matches_input_dir)

#     result = json_output.get_matches(match_data)

#     with open("tests/unit_tests/test_outputs/format_writer/matches/json/HAMAP.match.json", "w") as fh:
#         json.dump(result, fh)


# def test_mobidb_json_match(j_matches_input_dir):
#     match_data = load_match_data("MOBIDB", j_matches_input_dir)

#     result = json_output.get_matches(match_data)

#     with open("tests/unit_tests/test_outputs/format_writer/matches/json/MOBIDB.match.json", "w") as fh:
#         json.dump(result, fh)


# def test_panther_json_match(j_matches_input_dir):
#     match_data = load_match_data("PANTHER", j_matches_input_dir)

#     result = json_output.get_matches(match_data)

#     with open("tests/unit_tests/test_outputs/format_writer/matches/json/PANTHER.match.json", "w") as fh:
#         json.dump(result, fh)


# def test_phobius_json_match(j_matches_input_dir):
#     match_data = load_match_data("PHOBIUS", j_matches_input_dir)

#     result = json_output.get_matches(match_data)

#     with open("tests/unit_tests/test_outputs/format_writer/matches/json/PHOBIUS.match.json", "w") as fh:
#         json.dump(result, fh)


# def test_pirsf_json_match(j_matches_input_dir):
#     match_data = load_match_data("PIRSF", j_matches_input_dir)

#     result = json_output.get_matches(match_data)

#     with open("tests/unit_tests/test_outputs/format_writer/matches/json/PIRSF.match.json", "w") as fh:
#         json.dump(result, fh)


# def test_pirsr_json_match(j_matches_input_dir):
#     match_data = load_match_data("PIRSR", j_matches_input_dir)

#     result = json_output.get_matches(match_data)

#     with open("tests/unit_tests/test_outputs/format_writer/matches/json/PIRSR.match.json", "w") as fh:
#         json.dump(result, fh)


# def test_prints_json_match(j_matches_input_dir):
#     match_data = load_match_data("PRINTS", j_matches_input_dir)

#     result = json_output.get_matches(match_data)

#     with open("tests/unit_tests/test_outputs/format_writer/matches/json/PRINTS.match.json", "w") as fh:
#         json.dump(result, fh)


# def test_prosite_pattern_json_match(j_matches_input_dir):
#     match_data = load_match_data("PROSITE_PATTERNS", j_matches_input_dir)

#     result = json_output.get_matches(match_data)

#     with open("tests/unit_tests/test_outputs/format_writer/matches/json/PROSITE_PATTERNS.match.json", "w") as fh:
#         json.dump(result, fh)


# def test_prosite_profile_json_match(j_matches_input_dir):
#     match_data = load_match_data("PROSITE_PROFILES", j_matches_input_dir)

#     result = json_output.get_matches(match_data)

#     with open("tests/unit_tests/test_outputs/format_writer/matches/json/PROSITE_PROFILES.match.json", "w") as fh:
#         json.dump(result, fh)


# def test_sfld_json_match(j_matches_input_dir):
#     match_data = load_match_data("SFLD", j_matches_input_dir)

#     result = json_output.get_matches(match_data)

#     with open("tests/unit_tests/test_outputs/format_writer/matches/json/SFLD.match.json", "w") as fh:
#         json.dump(result, fh)


# def test_signalp_json_match(j_matches_input_dir):
#     match_data = load_match_data("SIGNALP", j_matches_input_dir)

#     result = json_output.get_matches(match_data)

#     with open("tests/unit_tests/test_outputs/format_writer/matches/json/SIGNALP.match.json", "w") as fh:
#         json.dump(result, fh)


# def test_smart_json_match(j_matches_input_dir):
#     match_data = load_match_data("SMART", j_matches_input_dir)

#     result = json_output.get_matches(match_data)

#     with open("tests/unit_tests/test_outputs/format_writer/matches/json/SMART.match.json", "w") as fh:
#         json.dump(result, fh)


# def test_superfamily_json_match(j_matches_input_dir):
#     match_data = load_match_data("SUPERFAMILY", j_matches_input_dir)

#     result = json_output.get_matches(match_data)

#     with open("tests/unit_tests/test_outputs/format_writer/matches/json/SUPERFAMILY.match.json", "w") as fh:
#         json.dump(result, fh)
