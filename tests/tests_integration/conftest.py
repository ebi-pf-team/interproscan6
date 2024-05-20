from pathlib import Path

import pytest


@pytest.fixture
def test_dir():
    return Path("tests/tests_integration/")


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
    return test_input_dir / "small_test.fasta"


@pytest.fixture
def expected_output_path(test_output_dir):
    return test_output_dir / "expected_output"


@pytest.fixture
def current_output_path(test_output_dir):
    return test_output_dir / "current_output"


@pytest.fixture
def applications():
    return "antifam,ncbifam,sfld,cdd,panther"


@pytest.fixture
def disable_precalc():
    return True
