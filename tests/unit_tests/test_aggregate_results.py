"""Test the python script that coordinates writing the output.

These test are intened to be run from the root of the repository using:
pytest -v
"""

import json

import pytest

from interproscan.scripts.output import aggregate_results


@pytest.fixture
def ar_input_paths(test_input_dir):
    input_dir = test_input_dir / "aggregate_results"
    cdd = input_dir / "cdd_matches.json"
    hamap = input_dir / "hamap_matches.json"
    return [cdd, hamap]


@pytest.fixture
def ar_expected_output(test_output_dir):
    _path =  test_output_dir / "aggregate_results/aggregated_results.json"
    with open(_path, "r") as fh:
        expected = json.load(fh)
    return expected


def test_aggregate_main(ar_input_paths, monkeypatch):
    def mock_aggregate_results(*args, **kwards):
        return

    input_paths = str([str(ar_input_paths[0]), str(ar_input_paths[1])])
    test_args = [
        "aggregate_results.py",
        input_paths,
        "output path"
    ]

    monkeypatch.setattr("sys.argv", test_args)
    monkeypatch.setattr(aggregate_results, "aggregate_results", mock_aggregate_results)

    assert aggregate_results.main() is None


def test_aggregating_results(ar_input_paths, ar_expected_output):
    assert ar_expected_output == aggregate_results.aggregate_results(ar_input_paths)
