"""Test the python script goterms.py

These test are intended to be run from the root of the repository using:
python -m pytest -v
"""

import json

import pytest

from interproscan.scripts.xrefs import goterms



@pytest.fixture
def goterms_path(test_input_dir):
    return test_input_dir / "xrefs/goterms"


@pytest.fixture
def matches_json_path(test_input_dir):
    return test_input_dir / "xrefs/match_entries.json"


@pytest.fixture
def expected_output(test_output_dir):
    _path = test_output_dir / "xref/entries_goterms.json"
    with open(_path, "r") as fh:
        output = json.load(fh)
    return output


@pytest.fixture
def temp_output(test_output_dir):
    dirpath = test_output_dir / "temp"
    filepath = dirpath / "unittest.ips6.json"
    dirpath.mkdir(parents=True, exist_ok=True)
    return filepath


def test_goterms_main(goterms_path, matches_json_path, temp_output, monkeypatch):
    def mock_add_goterms(*args, **kwards):
        return {}

    test_args = [
        "goterms.py",
        str(matches_json_path),
        str(goterms_path),
        str(temp_output)
    ]

    monkeypatch.setattr("sys.argv", test_args)
    monkeypatch.setattr(goterms, "add_goterms", mock_add_goterms)

    goterms.main()
    temp_output.unlink()


def test_add_goterms(goterms_path, expected_output, matches_json_path):
    def parse_goterms(goterm_list):
        goterms = {}
        for goterm_dict in goterm_list:
            go = goterm_dict["id"]
            if go not in goterms:
                goterms[go] = goterm_dict
            else:
                print(f"GO id {go} duplicated")
                goterms[go] = goterm_dict
        return goterms

    results = goterms.add_goterms(matches_json_path, str(goterms_path))

    for seq_id, seq_data in results.items():
        assert seq_id in expected_output, f"Seq id {seq_id} not in expected outputs"
        if seq_id in expected_output:
            for sig_acc, sig_data in seq_data.items():
                assert sig_acc in expected_output[seq_id], (
                    f"Signature acc {sig_acc} for seq id {seq_id} not in expected outputs"
                )
                if sig_acc in expected_output[seq_id]:
                    assert len(sig_data["entry"]["goXRefs"]) == len(expected_output[seq_id][sig_acc]["entry"]["goXRefs"]), (
                        f"Mismatch in the number of go terms for signature {sig_acc} for seq id {seq_id}"
                    )
                    result_goterms = parse_goterms(sig_data["entry"]["goXRefs"])
                    expected_goterms = parse_goterms(expected_output[seq_id][sig_acc]["entry"]["goXRefs"])
                    if len(result_goterms) > 0:
                        for go_id, go_data in result_goterms.items():
                            assert go_id in expected_goterms, (
                                f"GO term with ID {go_id} in result GO terms but not in "
                                f"the expected output for signature {sig_acc} for seq id {seq_id}"
                            )
                            if go_id in expected_goterms:
                                assert all(
                                    go_data[key] == expected_goterms[go_id][key]
                                    for key in ['name', 'databaseName', 'category']
                                ), (
                                    f"Mismatch in GO term data with GO id {go_id} "
                                    f"for signature {sig_acc} for seq id {seq_id}"
                                )
                        for go_id in expected_goterms:
                            assert go_id in result_goterms, (
                                f"GO term with ID {go_id} in expected GO terms but not in "
                                f"the result output for signature {sig_acc} for seq id {seq_id}"
                            )
