"""Test the python script that coordinates writing the output.

These test are intened to be run from the root of the repository using:
pytest -v
"""

import json

import pytest

from interproscan.scripts.output import aggregate_parsed_seqs


@pytest.fixture
def ar_input_paths(test_input_dir):
    input_dir = test_input_dir / "aggregate_parsed_seqs"
    path_1 = input_dir / "parse_seqs_1.json"
    path_2 = input_dir / "parse_seqs_2.json"
    return [path_1, path_2]


@pytest.fixture
def ar_expected_output(test_output_dir):
    _path =  test_output_dir / "aggregate_parsed_seqs/aggregate_parsed_seqs.json"
    with open(_path, "r") as fh:
        expected = json.load(fh)
    return expected


def test_aggregate_parse_seq_main(ar_input_paths, monkeypatch):
    def mock_aggregate(*args, **kwards):
        return {}

    input_paths = str([str(ar_input_paths[0]), str(ar_input_paths[1])])
    test_args = [
        "aggregate_parsed_seqs.py",
        input_paths,
        "test_output.json"
    ]

    monkeypatch.setattr("sys.argv", test_args)
    monkeypatch.setattr(aggregate_parsed_seqs, "aggregate_parsed_seqs", mock_aggregate)

    assert aggregate_parsed_seqs.main() is None


def test_aggregating_parsed_seqs(ar_input_paths, ar_expected_output):
    assert ar_expected_output == aggregate_parsed_seqs.aggregate_parsed_seqs(ar_input_paths)
