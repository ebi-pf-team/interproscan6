"""Test the python script that coordinates writing the output.

These test are intened to be run from the root of the repository using:
pytest -v
"""

import json

import pytest

from interproscan.scripts.output import aggregate_results


@pytest.fixture
def ar_input_path(test_input_dir):
    return test_input_dir / "aggregate_results/cdd_matches.json"


@pytest.fixture
def ar_expected_output_path(test_output_dir):
    return test_output_dir / "aggregate_results/aggregated_results.json"


@pytest.fixture
def temp_path(test_output_dir):
    dirpath = test_output_dir / "temp"
    filepath = dirpath / "temp.aggregated.results.json"
    dirpath.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w") as fh:
        json.dump({}, fh)
    return filepath


def test_aggregate_main(ar_input_path, temp_path, monkeypatch):
    def mock_aggregate_results(*args, **kwards):
        return

    test_args = [
        "aggregate_results.py",
        str(temp_path),
        str(ar_input_path),
    ]

    monkeypatch.setattr("sys.argv", test_args)
    monkeypatch.setattr(aggregate_results, "aggregate_results", mock_aggregate_results)

    assert aggregate_results.main() is None

    temp_path.unlink()


def test_aggregating_results(ar_input_path, ar_expected_output_path):
    with open(ar_expected_output_path, "r") as fh:
        expected = json.load(fh)

    assert expected == aggregate_results.aggregate_results(
        str(ar_expected_output_path),
        str(ar_input_path)
    )
