"""Test the python script that generates the MD5 hash dict of seqs.

These test are intened to be run from the root of the repository using:
pytest -v
"""

import json

import pytest

from interproscan.scripts.parse_sequence import parse_sequence


@pytest.fixture
def temp_path(test_output_dir):
    dirpath = test_output_dir / "temp"
    filepath = dirpath / "temp.parsed.seqs.json"
    dirpath.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w") as fh:
        json.dump({}, fh)
    return filepath


def test_parse_seq_main(temp_path, monkeypatch):
    def mock_parse_seqs(*args, **kwards):
        return

    test_args = [
        "parse_sequence.py",
        "fake_path_to_FASTA",
        "fake_path_to_FASTA",
        "false",
        "antifam,sfld",
        str(temp_path)
    ]

    monkeypatch.setattr("sys.argv", test_args)
    monkeypatch.setattr(parse_sequence, "parse", mock_parse_seqs)
    
    assert parse_sequence.main() is None
