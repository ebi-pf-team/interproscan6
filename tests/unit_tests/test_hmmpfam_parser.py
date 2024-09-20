"""Test the python script that parses the HMMER pfam output, parse_hmmpfam_out.py

These test are intended to be run from the root of the repository using:
python -m pytest -v
"""

import pytest

from interproscan.scripts.hmmer import parse_hmmpfam_out


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
