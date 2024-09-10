"""Test the python script paint_annotations.py

These test are intended to be run from the root of the repository using:
python -m pytest -v
"""

import json

import pytest

from interproscan.scripts.xrefs import paint_annotations


@pytest.fixture
def paint_annotations_path(test_input_dir):
    
    return test_input_dir / "xrefs/paint_annotations"


@pytest.fixture
def matches_json_path(test_input_dir):
    return test_input_dir / "xrefs/paint_annotations/entries_matches.json"


@pytest.fixture
def expected_output(test_output_dir):
    _path = test_output_dir / "xref/paint_annotations.json"
    with open(_path, "r") as fh:
        output = json.load(fh)
    return output


@pytest.fixture
def temp_output(test_output_dir):
    dirpath = test_output_dir / "temp"
    filepath = dirpath / "unittest.ips6.json"
    dirpath.mkdir(parents=True, exist_ok=True)
    return filepath


def test_paint_annotations_main(paint_annotations_path, matches_json_path, temp_output, monkeypatch):
    def add_paint_annotations(*args, **kwards):
        return {}

    test_args = [
        "paint_annotations.py",
        str(matches_json_path),
        str(paint_annotations_path),
        str(temp_output)
    ]

    monkeypatch.setattr("sys.argv", test_args)
    monkeypatch.setattr(paint_annotations, "add_paint_annotations", add_paint_annotations)

    paint_annotations.main()
    temp_output.unlink()


def test_add_paint_annotations(paint_annotations_path, expected_output, matches_json_path):
    def parse_paint_annotations(goterm_list):
        paint_annotations = {}
        for pathway_dict in goterm_list:
            pw = pathway_dict["id"]
            if pw not in paint_annotations:
                paint_annotations[pw] = pathway_dict
            else:
                print(f"Pathway id {pw} duplicated")
                paint_annotations[pw] = pathway_dict
        return paint_annotations

    results = paint_annotations.add_paint_annotations(matches_json_path, paint_annotations_path)

    for seq_id, seq_data in results.items():
        assert seq_id in expected_output, f"Seq id {seq_id} not in expected outputs"
        if seq_id in expected_output:
            for sig_acc, sig_data in seq_data.items():
                assert sig_acc in expected_output[seq_id], (
                    f"Signature acc {sig_acc} for seq id {seq_id} not in expected outputs"
                )
                if sig_acc in expected_output[seq_id]:
                    assert sig_data["proteinClass"] == expected_output[seq_id][sig_acc]["proteinClass"], (
                        f"Mismatched proteinClass for signature {sig_acc} for seq id {seq_id}"
                    )
