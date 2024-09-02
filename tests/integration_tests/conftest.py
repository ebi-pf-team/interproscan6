from pathlib import Path

import pytest


@pytest.fixture
def test_dir():
    return Path('tests/integration_tests/')


@pytest.fixture
def test_input_dir(test_dir):
    dir_path = test_dir / 'test_inputs'
    return dir_path


@pytest.fixture
def test_output_dir(test_dir):
    dir_path = test_dir / 'test_outputs'
    return dir_path


@pytest.fixture
def input_path(test_input_dir):
    return test_input_dir / 'integration_input.fasta'


@pytest.fixture
def output_path(test_output_dir):
    return test_output_dir / 'integration_input.fasta.ips6'


@pytest.fixture
def expected_output_path(test_output_dir):
    return test_output_dir / 'expected_output_representative'


@pytest.fixture
def applications():
    return 'CDD,NCBIfam,Pfam,PROSITE_Profiles,SMART'


@pytest.fixture
def disable_precalc():
    return True
