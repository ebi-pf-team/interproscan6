"""Test the python script entries.py

These test are intended to be run from the root of the repository using:
python -m pytest -v
"""

import json

import pytest

from interproscan.scripts.xrefs import entries


@pytest.fixture
def entries_path(test_input_dir):
    return test_input_dir / "xrefs/entries.json"


@pytest.fixture
def matches_json_path(test_input_dir):
    return test_input_dir / "xrefs/subworkflow/xref_matches_input.json"


@pytest.fixture
def expected_output(test_output_dir):
    _path = test_output_dir / "xref/match_entries.json"
    with open(_path, "r") as fh:
        output = json.load(fh)
    return output


@pytest.fixture
def temp_output(test_output_dir):
    dirpath = test_output_dir / "temp"
    filepath = dirpath / "unittest.json"
    dirpath.mkdir(parents=True, exist_ok=True)
    return filepath


def test_entries_main(entries_path, matches_json_path, temp_output, monkeypatch):
    def mock_add_entries(*args, **kwards):
        return {}

    test_args = [
        "entries.py",
        str(matches_json_path),
        str(entries_path),
        str(temp_output)
    ]

    monkeypatch.setattr("sys.argv", test_args)
    monkeypatch.setattr(entries, "add_entries", mock_add_entries)

    entries.main()
    temp_output.unlink()


def test_add_entries(entries_path, expected_output, matches_json_path):
    results = entries.add_entries(matches_json_path, entries_path)

    for seq_id, seq_data in results.items():
        assert seq_id in expected_output, f"Seq id {seq_id} not in expected outputs"
        if seq_id in expected_output:
            for sig_acc, sig_data in seq_data.items():
                assert sig_acc in expected_output[seq_id], (
                    f"Signature acc {sig_acc} for seq id {seq_id} not in expected outputs"
                )
                if sig_acc in expected_output[seq_id]:
                    assert all(
                        sig_data["entry"][key] == expected_output[seq_id][sig_acc]["entry"][key]
                        for key in [
                            'accession', 'name', 'description', 'type',
                            'goXRefs', 'pathwayXRefs'
                        ]
                    ), (
                        f"Mismatch in entries data for signature acc {sig_acc} "
                        f"for seq id {seq_id}"
                    )
