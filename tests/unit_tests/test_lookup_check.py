"""Test the python script scripts/lookup/lookup_check.py

These test are intened to be run from the root of the repository using:
python -m pytest -v
"""

import json
import pytest
import urllib.request

from interproscan.scripts.lookup import lookup_check


@pytest.fixture
def lookup_check_input_dir(test_input_dir):
    return test_input_dir / "precalc_match_lookup"


@pytest.fixture
def lookup_check_out_dir(test_output_dir):
    return test_output_dir / "precalc_match_lookup"


@pytest.fixture
def lookup_outfile_path(test_output_dir):
    return test_output_dir / "temp/lookup.out.json"


@pytest.fixture
def parsed_seqs_path(lookup_check_input_dir):
    return lookup_check_input_dir / "parsed_sequences"


def test_lookup_check_main(
    parsed_seqs_path, lookup_check_out_dir,
    lookup_outfile_path, capsys, monkeypatch
):
    def mock_check_precalc(*args, **kwards):
        return (["MD5_1", "MD5_2"], None)

    test_args = [
        "lookup_check",
        str(parsed_seqs_path),
        "fake_url",
        3,
        lookup_outfile_path
    ]

    monkeypatch.setattr("sys.argv", test_args)
    monkeypatch.setattr(lookup_check, "check_precalc", mock_check_precalc)

    lookup_check.main()

    lookup_outfile_path.unlink()


def test_check_precalc(monkeypatch):
    class MockResponse:
        def read(self):
            return b"MD5_1\nMD5_2\nMD5_3"

    def mock_urlopen(url):
        return MockResponse()

    monkeypatch.setattr(urllib.request, "urlopen", mock_urlopen)

    md5_list = ["MD5_1", "MD5_2", "MD5_3", "MD5_4"]
    url = "https://fake_url"
    result = lookup_check.check_precalc(md5_list, url, retries=3)

    assert result == (["MD5_1", "MD5_2", "MD5_3"], None)
