"""Test the python script that tests if a FASTA file contains any non-nucleic sequences

These test are intened to be run from the root of the repository using:
pytest -v
"""

import pytest

from scripts.pre_checks import check_nucleic_seq


@pytest.fixture
def precheck_input_dir(test_input_dir):
    return test_input_dir / "pre_analysis_checks"


def test_is_nucleic_true():
    """Test is_nucleic when a nucleic seq is provided"""
    seq = "ATGC"
    assert check_nucleic_seq.is_nucleic(seq)


def test_is_nucleic_false():
    """Test is_nucleic when a NON-nucleic seq is provided"""
    seq = "ATGCPROTEIN"
    assert check_nucleic_seq.is_nucleic(seq) is False


def test_provided_nucleic(precheck_input_dir):
    """Test main() when a nucleic seq is provided"""
    _path = precheck_input_dir / "nt_seqs.fasta"
    check_nucleic_seq.main(_path)


def test_provided_non_nucleic(precheck_input_dir):
    """Test main() when a NON-nucleic seq is provided"""
    _path = precheck_input_dir / "protein_seqs.fasta"
    with pytest.raises(check_nucleic_seq.NuleicError) as pytest_wrapped_e:
        check_nucleic_seq.main(_path)
    assert pytest_wrapped_e.type == check_nucleic_seq.NuleicError
