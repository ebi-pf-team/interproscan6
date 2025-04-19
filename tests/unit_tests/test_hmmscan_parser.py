"""Test the python script that parses the HMMER pfam output, parser_scan_out.py

These test are intended to be run from the root of the repository using:
python -m pytest -v
"""

import json

import pytest

from interproscan.scripts.hmmer import parser_scan_out


@pytest.fixture
def hmmscan_out(test_input_dir):
    return test_input_dir / "hmmscan_parser/hmmscan.out"


@pytest.fixture
def expected_parser_output(test_output_dir):
    _path = test_output_dir / "hmmscan_parser/expected_parser_output.json"
    with open(_path, "r") as fh:
        data = json.load(fh)
    return data


@pytest.fixture
def expected_matches_output(test_output_dir):
    _path = test_output_dir / "hmmscan_parser/expected_matches_output.json"
    with open(_path, "r") as fh:
        data = json.load(fh)
    return data


@pytest.fixture
def temp_output(test_output_dir):
    dirpath = test_output_dir / "temp"
    filepath = dirpath / "unittest.json"
    dirpath.mkdir(parents=True, exist_ok=True)
    return filepath


def test_hmmscan_main(monkeypatch, temp_output):
    def mock_parser_scan_out(*args, **kwards):
        return {}

    test_args = [
        "parser_scan_out.py",
        "hmmscan_outfile",
        "member_db",
        temp_output
    ]

    monkeypatch.setattr("sys.argv", test_args)
    monkeypatch.setattr(parser_scan_out, "parse", mock_parser_scan_out)

    parser_scan_out.main()

    temp_output.unlink()


def test_parsing_hmmer2(expected_matches_output, hmmscan_out, monkeypatch):
    def mock_add_match(*args, **kwards):
        return expected_matches_output

    monkeypatch.setattr(parser_scan_out, "add_match", mock_add_match)

    member_db = "unit-test"

    assert expected_matches_output == parser_scan_out.parse(
        hmmscan_out,
        member_db
    )


def test_add_hmmscan_matches(expected_matches_output):
    protein = parser_scan_out.QueryProtein()
    query_line = "Query:       sp|A2SLW2|1A1D_METPP  [L=338]"
    query_match = parser_scan_out.ACC_LINE.match(query_line.strip())
    protein.get_seq_id(query_match)
    # model 1
    current_model = "PIRSF006278"
    sig_match = "   7.2e-127  420.8   0.0   7.9e-127  420.7   0.0    1.0  1  PIRSF006278"
    model_line_pattern = parser_scan_out.MODEL_HIT_LINE.match(sig_match.strip())
    protein.get_model_data(model_line_pattern)
    domain_match = "   1 !  420.7   0.0  9.7e-130  7.9e-127      10     349 ..       2     335 ..       1     338 [] 0.98"
    domain_line_pattern = parser_scan_out.DOMAIN_HIT_LINE.match(domain_match.strip())
    protein.get_domain_data(current_model, domain_line_pattern)

    # model 2
    current_model = "PIRSF500824"
    sig_match = "    0.00012   18.4   0.0       0.14    8.2   0.0    2.3  3  PIRSF500824"
    model_line_pattern = parser_scan_out.MODEL_HIT_LINE.match(sig_match.strip())
    protein.get_model_data(model_line_pattern)
    domain_match = "   1 !    8.2   0.0   0.00017      0.14      79     157 ..      14      96 ..       8     128 .. 0.65"
    domain_line_pattern = parser_scan_out.DOMAIN_HIT_LINE.match(domain_match.strip())
    protein.get_domain_data(current_model, domain_line_pattern)
    domain_match = "   2 !    7.8   0.0   0.00024      0.19     245     305 ..     178     234 ..     172     253 .. 0.87"
    domain_line_pattern = parser_scan_out.DOMAIN_HIT_LINE.match(domain_match.strip())
    protein.get_domain_data(current_model, domain_line_pattern)
    domain_match = "   3 ?   -4.0   0.0      0.88   7.2e+02     365     384 ..     273     292 ..     270     311 .. 0.75"
    domain_line_pattern = parser_scan_out.DOMAIN_HIT_LINE.match(domain_match.strip())
    protein.get_domain_data(current_model, domain_line_pattern)

    matches = {}
    member_db = "unit-test"

    assert expected_matches_output == parser_scan_out.add_match(matches, protein, member_db)
