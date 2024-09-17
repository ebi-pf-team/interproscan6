"""Test the python script that generates the MD5 hash dict of seqs.

These test are intened to be run from the root of the repository using:
pytest -v
"""

import json

import pytest

from interproscan.scripts.parse_sequence import parse_sequence


@pytest.fixture
def temp_path(test_output_dir):
    dirpath = test_output_dir / "temp"
    filepath = dirpath / "temp.parsed.seqs.json"
    dirpath.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w") as fh:
        json.dump({}, fh)
    return filepath


@pytest.fixture
def input_illegal_fasta_path(test_input_dir):
    return test_input_dir / "parse_sequences/illegal.faa"


@pytest.fixture
def input_legal_fasta_path(test_input_dir):
    return test_input_dir / "parse_sequences/allowed.faa"


def test_parse_seq_main(temp_path, monkeypatch):
    def mock_parse_seqs(*args, **kwards):
        return

    test_args = [
        "parse_sequence.py",
        "fake_path_to_FASTA",
        "fake_path_to_FASTA",
        "false",
        "antifam,sfld",
        str(temp_path)
    ]

    monkeypatch.setattr("sys.argv", test_args)
    monkeypatch.setattr(parse_sequence, "parse", mock_parse_seqs)

    assert parse_sequence.main() is None


def test_store_nucleic_seq():
    seq_obj = parse_sequence.Sequence()
    passing_nucleic = True
    seq_obj.get_seq_key(">Bob", passing_nucleic)
    seq_obj.get_seq("atgaaatttccccaaaggggaaa")
    sequences = {}
    expected = {"Bob": seq_obj}
    assert expected == parse_sequence.store_seq(
        seq_obj, sequences, passing_nucleic=passing_nucleic
    )


def test_store_protein_seq():
    passing_nucleic = True
    nucleic_seq_obj = parse_sequence.Sequence()
    nucleic_seq_obj.get_seq_key(">Bob", passing_nucleic)
    nucleic_seq_obj.get_seq("atgaaatttccccaaaggggaaa")
    seq_obj = parse_sequence.Sequence()
    passing_nucleic = False
    seq_obj.get_seq_key(">orf1 source=Bob coords=143..286 length=48 frame=2 desc=", passing_nucleic)
    seq_obj.get_seq("MAMAMAMAMAMAMAMAMAMAMLAMA")
    nucleic_seqs = {
        "Bob": nucleic_seq_obj
    }
    sequences = {}
    expected = {
        'orf1': {
            'seq_id': 'orf1 source=Bob coords=143..286 length=48 frame=2 desc=',
            'name': None,
            'sequence': 'MAMAMAMAMAMAMAMAMAMAMLAMA',
            'md5': '95ccc5915e0f68f7f43cfebda07fb866',
            'length': 25,
            'nt_seq_id': 'Bob',
            'nt_name': 'Bob',
            'nt_sequence': 'atgaaatttccccaaaggggaaa',
            'nt_md5': '29e4b52298105869362604ee4938bb92'
        }
    }
    assert expected == parse_sequence.store_seq(
        seq_obj, sequences, nucleic_seqs, passing_nucleic=passing_nucleic
    )


def test_parse_no_illegal_chars(input_legal_fasta_path, monkeypatch):
    def mock_store_seq(*args, **kwards):
        seq_obj = parse_sequence.Sequence()
        passing_nucleic = False
        seq_obj.get_seq_key(">TestProteinId", passing_nucleic)
        seq_obj.get_seq("MAMAMAMAMAMAMAMAMAMAMLAMA")
        return {"TestProteinId": seq_obj}

    seq_obj = parse_sequence.Sequence()
    passing_nucleic = False
    seq_obj.get_seq_key(">TestProteinId", passing_nucleic)
    seq_obj.get_seq("MAMAMAMAMAMAMAMAMAMAMLAMA")
    monkeypatch.setattr(parse_sequence, "store_seq", mock_store_seq)

    parse_sequence.parse(
        input_legal_fasta_path,
        "antifam,sfld",
    )


def test_parse_illegal_chars(input_illegal_fasta_path, monkeypatch):
    def mock_store_seq(*args, **kwargs):
        seq_obj = parse_sequence.Sequence()
        passing_nucleic = False
        seq_obj.get_seq_key(">TestProteinId", passing_nucleic)
        seq_obj.get_seq("MAMAMAMA*.__<>MAMAMAMAMAMAMLAMA")
        return {"TestProteinId": seq_obj}

    seq_obj = parse_sequence.Sequence()
    passing_nucleic = False
    seq_obj.get_seq_key(">TestProteinId", passing_nucleic)
    seq_obj.get_seq("MAMAMAMA*.__<>MAMA$£%^AMAMAMAMLAMAA")
    monkeypatch.setattr(parse_sequence, "store_seq", mock_store_seq)

    with pytest.raises(parse_sequence.IllegalCharError) as pytest_wrapped_e:
        parse_sequence.parse(
            input_illegal_fasta_path,
            "antifam,sfld",
        )
    assert pytest_wrapped_e.type == parse_sequence.IllegalCharError
    assert "Illegal characters detected" in str(pytest_wrapped_e.value)
    assert "Sequence TESTESEMLOOKUP contains illegal character(s)" in str(pytest_wrapped_e.value)
