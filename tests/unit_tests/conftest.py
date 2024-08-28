from pathlib import Path

import sys

import pytest


format_writer_path = Path(__file__).resolve().parent.parent.parent / 'interproscan' / 'scripts' / 'format_writer'
lookup_path = Path(__file__).resolve().parent.parent.parent / 'interproscan' / 'scripts' / 'lookup'
output_path = Path(__file__).resolve().parent.parent.parent / 'interproscan' / 'scripts' / 'output'

sys.path.insert(0, str(format_writer_path))
sys.path.insert(0, str(lookup_path))
sys.path.insert(0, str(output_path))


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
