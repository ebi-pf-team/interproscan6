from pathlib import Path

import pytest


@pytest.fixture
def test_input_dir():
    dir_path = Path('tests/data')
    return dir_path


@pytest.fixture
def test_output_dir():
    dir_path = Path('tests/integration_tests')
    return dir_path


@pytest.fixture
def input_path(test_input_dir):
    return test_input_dir / 'test_prot.fa'


@pytest.fixture
def expected_result_path(test_input_dir):
    return test_input_dir / 'iprscan5.zip'


@pytest.fixture
def output_path(test_output_dir):
    return test_output_dir / 'test_prot.fa.ips6'


@pytest.fixture
def data_dir():
    return Path('data')


@pytest.fixture
def disable_precalc():
    return True


def pytest_addoption(parser):
    parser.addoption(
        "--applications",
        action="store",
        default="AntiFam,CDD,Coils,FunFam,Gene3D,HAMAP,MobiDB_lite,NCBIfam,Panther,Pfam,PIRSF,PIRSR,PRINTS,PROSITE_PATTERNS,PROSITE_Profiles,SFLD,SMART,SUPERFAMILY",
    )


@pytest.fixture
def applications(request):
    return request.config.getoption("--applications")
