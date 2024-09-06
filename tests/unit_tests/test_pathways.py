"""Test the python script pathways.py

These test are intended to be run from the root of the repository using:
python -m pytest -v
"""

import json

import pytest

from interproscan.scripts.xrefs import pathways


@pytest.fixture
def pathways_path(test_input_dir):
    return test_input_dir / "xrefs/goterms_pathways/pathways"


@pytest.fixture
def matches_json_path(test_input_dir):
    return test_input_dir / "xrefs/goterms_pathways/match_entries.json"


@pytest.fixture
def expected_output(test_output_dir):
    _path = test_output_dir / "xref/entries_pathways.json"
    with open(_path, "r") as fh:
        output = json.load(fh)
    return output


@pytest.fixture
def temp_output(test_output_dir):
    dirpath = test_output_dir / "temp"
    filepath = dirpath / "unittest.ips6.json"
    dirpath.mkdir(parents=True, exist_ok=True)
    return filepath


def test_pathways_main(pathways_path, matches_json_path, temp_output, monkeypatch):
    def mock_add_pathways(*args, **kwards):
        return {}

    test_args = [
        "pathways.py",
        str(matches_json_path),
        str(pathways_path),
        str(temp_output)
    ]

    monkeypatch.setattr("sys.argv", test_args)
    monkeypatch.setattr(pathways, "add_pathways", mock_add_pathways)

    pathways.main()
    temp_output.unlink()


def test_add_pathways(pathways_path, expected_output, matches_json_path):
    def parse_pathways(goterm_list):
        pathways = {}
        for pathway_dict in goterm_list:
            pw = pathway_dict["id"]
            if pw not in pathways:
                pathways[pw] = pathway_dict
            else:
                print(f"Pathway id {pw} duplicated")
                pathways[pw] = pathway_dict
        return pathways

    results = pathways.add_pathways(matches_json_path, str(pathways_path))

    for seq_id, seq_data in results.items():
        assert seq_id in expected_output, f"Seq id {seq_id} not in expected outputs"
        if seq_id in expected_output:
            for sig_acc, sig_data in seq_data.items():
                assert sig_acc in expected_output[seq_id], (
                    f"Signature acc {sig_acc} for seq id {seq_id} not in expected outputs"
                )
                if sig_acc in expected_output[seq_id]:
                    assert len(sig_data["entry"]["pathwayXRefs"]) == len(expected_output[seq_id][sig_acc]["entry"]["pathwayXRefs"]), (
                        f"Mismatch in the number of pathways for signature {sig_acc} for seq id {seq_id}"
                    )
                    result_goterms = parse_pathways(sig_data["entry"]["pathwayXRefs"])
                    expected_pathways = parse_pathways(expected_output[seq_id][sig_acc]["entry"]["pathwayXRefs"])
                    if len(result_goterms) > 0:
                        for pw_id, pw_data in result_goterms.items():
                            assert pw_id in expected_pathways, (
                                f"Pathway with ID {pw_id} in result Pathways but not in "
                                f"the expected output for signature {sig_acc} for seq id {seq_id}"
                            )
                            if pw_id in expected_pathways:
                                assert all(
                                    pw_data[key] == expected_pathways[pw_id][key]
                                    for key in ['name', 'databaseName']
                                ), (
                                    f"Mismatch in Pathway data with Pathway id {pw_id} "
                                    f"for signature {sig_acc} for seq id {seq_id}"
                                )
                        for pw_id in expected_pathways:
                            assert pw_id in result_goterms, (
                                f"Pathway with ID {pw_id} in expected Pathways but not in "
                                f"the result output for signature {sig_acc} for seq id {seq_id}"
                            )
