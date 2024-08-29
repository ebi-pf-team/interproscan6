"""Test the python script that coordinates writing the output.

These test are intened to be run from the root of the repository using:
python -m pytest -v
"""

import pytest

from interproscan.scripts.output import write_results


@pytest.fixture
def wr_input_dir(test_input_dir):
    return test_input_dir / "write_results"


def test_wr_main(monkeypatch):
    def mock_write_results(*args, **kwards):
        return

    test_args = [
        "write_results.py",
        "sequences",
        "matches",
        "formats,string",
        "outfile_name",
        "version",
        "true"
    ]

    monkeypatch.setattr("sys.argv", test_args)
    monkeypatch.setattr(write_results, "write_results", mock_write_results)

    assert write_results.main() is None


def test_wr_write_results_nucleic(wr_input_dir, monkeypatch):
    """Test write_results() when the input is nucleic sequences"""
    def mock_write_results(*args, **kwards):
        return

    monkeypatch.setattr(write_results, "tsv_output", mock_write_results)
    monkeypatch.setattr(write_results, "tsv_pro_output", mock_write_results)
    monkeypatch.setattr(write_results, "build_json_output_nucleic", mock_write_results)
    monkeypatch.setattr(write_results, "build_xml_output_nucleic", mock_write_results)

    matches_json = wr_input_dir / "nucleic/matches.json"
    parsed_seq_json = wr_input_dir / "nucleic/parsed_seqs.json"
    out_format = ["TSV", "JSON", "XML", "TSV-PRO"]
    filename = "test_write_results.fasta"
    version = "unit.test"
    nucleic = True

    assert write_results.write_results(
        parsed_seq_json,
        matches_json,
        out_format,
        filename,
        version,
        nucleic
    ) is None


def test_wr_write_results_protein(wr_input_dir, monkeypatch):
    """Test write_results() when the input is protein sequences"""
    def mock_write_results(*args, **kwards):
        return

    monkeypatch.setattr(write_results, "tsv_output", mock_write_results)
    monkeypatch.setattr(write_results, "tsv_pro_output", mock_write_results)
    monkeypatch.setattr(write_results, "build_json_output_protein", mock_write_results)
    monkeypatch.setattr(write_results, "build_xml_output_protein", mock_write_results)

    matches_json = wr_input_dir / "protein/matches.json"
    parsed_seq_json = wr_input_dir / "protein/parsed_seqs.json"
    out_format = ["TSV", "JSON", "XML", "TSV-PRO"]
    filename = "test_write_results.fasta"
    version = "unit.test"
    nucleic = False

    assert write_results.write_results(
        parsed_seq_json,
        matches_json,
        out_format,
        filename,
        version,
        nucleic
    ) is None
