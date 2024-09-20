"""Test the python script that parses the HMMER pfam output, parse_hmmpfam_out.py

These test are intended to be run from the root of the repository using:
python -m pytest -v
"""

import json

import pytest

from interproscan.scripts.hmmer import parse_hmmpfam_out


@pytest.fixture
def expected_matches_output(test_input_dir):
    _path = test_input_dir / "hmmer2_parser/expected_matches_output.json"
    with open(_path, "r") as fh:
        data = json.load(fh)
    return data


@pytest.fixture
def temp_output(test_output_dir):
    dirpath = test_output_dir / "temp"
    filepath = dirpath / "unittest.ips6.json"
    dirpath.mkdir(parents=True, exist_ok=True)
    return filepath


def test_hmmpfam_main(monkeypatch, temp_output):
    def mock_parse_hmmpfam_out(*args, **kwards):
        return {}

    test_args = [
        "parse_hmmpfam_out.py",
        "hmmpfam_outfile",
        "member_db",
        temp_output
    ]

    monkeypatch.setattr("sys.argv", test_args)
    monkeypatch.setattr(parse_hmmpfam_out, "parse_hmmpfam_out", mock_parse_hmmpfam_out)

    parse_hmmpfam_out.main()

    temp_output.unlink()


def test_add_hmmpfam_matches(expected_matches_output):
    protein = parse_hmmpfam_out.QueryProtein()
    query_line = "Query sequence: UPI0010EB9177"
    protein.get_seq_id(query_line)
    sig_match = "SM00870                                                 462.1   2.8e-134   1"
    model_line_pattern = parse_hmmpfam_out.SIGNATURE_DATA_LINE.match(sig_match)
    protein.get_model_data(model_line_pattern)
    domain_match = "SM00870          1/1      51   364 ..     1   408 []   462.1 2.8e-134"
    domain_line_pattern = parse_hmmpfam_out.DOMAIN_DATA_LINE.match(domain_match)
    protein.get_domain_data(domain_line_pattern)

    matches = {}
    member_db = "unit-test"

    assert expected_matches_output == parse_hmmpfam_out.add_match(matches, protein, member_db)
