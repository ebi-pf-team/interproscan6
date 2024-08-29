"""Test the python script that coordinates writing the output.

These test are intened to be run from the root of the repository using:
python -m pytest -v
"""

import json

import pytest

from interproscan.scripts.output.format_writer import json_output


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


@pytest.fixture
def j_nucleic_seq_matches(test_input_dir):
    _path = test_input_dir / "format_writer/nucleic/seq_match_dict.json"
    with open(_path, "r") as fh:
        _input = json.load(fh)
    return _input


@pytest.fixture
def j_expected_protein(test_output_dir):
    _path = test_output_dir / "format_writer/build_funcs/expected_json_protein.json"
    with open(_path, "r") as fh:
        expected_output = json.load(fh)
    return expected_output


@pytest.fixture
def j_expected_nucleic(test_output_dir):
    _path = test_output_dir / "format_writer/build_funcs/expected_json_nucleic.json"
    with open(_path, "r") as fh:
        expected_output = json.load(fh)
    return expected_output


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
