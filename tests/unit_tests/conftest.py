from pathlib import Path

import sys

import pytest


sys.path.insert(
    0,
    Path(__file__).resolve().parent.parent.parent / 'interproscan' / 'scripts'
)


@pytest.fixture
def test_dir():
    return Path("tests/unit_tests")


@pytest.fixture
def test_input_dir(test_dir):
    dir_path = test_dir / "test_inputs"
    return dir_path


@pytest.fixture
def test_output_dir(test_dir):
    dir_path = test_dir / "test_outputs"
    return dir_path
