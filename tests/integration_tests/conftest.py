from pathlib import Path

import pytest


@pytest.fixture
def test_dir():
    return Path("tests/integration_tests/")


@pytest.fixture
def test_input_dir(test_dir):
    dir_path = test_dir / "test_inputs"
    return dir_path


@pytest.fixture
def test_output_dir(test_dir):
    dir_path = test_dir / "test_outputs"
    return dir_path


@pytest.fixture
def input_path(test_input_dir):
    return test_input_dir / "pfam_fragm.fasta"


@pytest.fixture
def current_output_path(test_output_dir):
    return test_output_dir / "current_output"


@pytest.fixture
def expected_output_path(test_output_dir):
    return test_output_dir / "expected_output_prosite_patterns"


@pytest.fixture
def applications():
    return "prosite_patterns"


@pytest.fixture
def disable_precalc():
    return True
