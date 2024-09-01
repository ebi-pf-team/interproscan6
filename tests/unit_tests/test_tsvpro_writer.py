"""Test the python script that coordinates writing the TSV-PRO output.

These test are intened to be run from the root of the repository using:
python -m pytest -v
"""

import json

from pathlib import Path

import pytest

from interproscan.scripts.output.format_writer import tsv_output



def load_tsvpro_results(tsv_path: Path) -> dict[str, dict[str, str]]:
    results = {}
    with open(tsv_path, "r") as fh:
        for line in fh:
            data = line.split()
            seq_id = data[3]
            mem_db = data[0]
            sig_acc = data[4]
            start, end = data[6], data[7]
            release = data[1]
            fragments = data[8]
            score = data[10]
            hmm_start, hmm_end = data[11], data[12]
            hmm_len = data[13]
            location_score = data[16]
            location_evalue = data[17]
            if mem_db not in results:
                results[mem_db] = {}
            if seq_id not in results[mem_db]:
                results[seq_id][mem_db] = {}
            if mem_db not in results[seq_id][mem_db]:
                results[seq_id][mem_db][sig_acc] = {}
            position = f"{start}-{end}--{fragments}"
            if position not in results[seq_id][mem_db][sig_acc]:
                results[seq_id][mem_db][sig_acc][position] = {
                    'release': release,
                    'fragments': fragments,
                    'score': score,
                    'hmm_start': hmm_start,
                    'hmm_end': hmm_end,
                    'hmm_len': hmm_len,
                    'location_score': location_score,
                    'location_evalue': location_evalue
                }
    return results


@pytest.fixture
def tsvout_path(test_output_dir):
    dirpath = test_output_dir / "temp"
    filepath = dirpath / "tsv.unittest.ips6.tsv"
    dirpath.mkdir(parents=True, exist_ok=True)
    return filepath


@pytest.fixture
def expected_tsvpro_outdir(test_output_dir):
    return test_output_dir / "format_writer/tsv-pro"


@pytest.fixture
def tsv_seq_match_dir(test_input_dir):
    return test_input_dir / "format_writer/tsv"


def load_seq_matches_dict(member_db, seq_match_dir):
    with open((seq_match_dir / f"{member_db}.seq-match-dict.json"), "r") as fh:
        seq_match_dict = json.load(fh)
    return seq_match_dict


def compare_tsvpro_output(
    test_tsv_out: dict[str, dict[str, dict[str, str]]],
    expected_tsvpro_out:  dict[str, dict[str, dict[str, str]]],
) -> None:
    for seq_id, member_dbs in test_tsv_out.items():
        assert seq_id in expected_tsvpro_out, f"seq_id {seq_id} not in expected_tsvpro_out"
        for member_db, sig_accs in member_dbs.items():
            assert member_db in expected_tsvpro_out[seq_id], \
                f"member_db {member_db} not in expected_tsvpro_out[{seq_id}]"
            for sig_acc, positions in sig_accs.items():
                assert sig_acc in expected_tsvpro_out[seq_id][member_db], \
                    f"sig_acc {sig_acc} not in expected_tsvpro_out[{seq_id}][{member_db}]"
                for position, match in positions.items():
                    expected_match = expected_tsvpro_out[seq_id][member_db][sig_acc][position]
                    assert position in expected_tsvpro_out[seq_id][member_db][sig_acc], \
                        f"position {position} not in expected_tsvpro_out[{seq_id}][{member_db}][{sig_acc}]"
                    assert all(
                        match[key] == expected_match[key]
                        for key in [
                            'release', 'fragments', 'score', 'hmm_start', 'hmm_end',
                            'hmm_len', 'location_score', 'location_evalue'
                        ]
                    ), f"Mismatch in match details for {seq_id}, {member_db}, {sig_acc}, {position}"



def test_tsvpro_antifam_output(tsvout_input, tsvout_path, expected_tsvpro_outdir):
    member_db = "ANTIFAM"
    tsv_output.tsv_pro_output(tsvout_input, (expected_tsvpro_outdir / f"{member_db}.tsvpro.tsv"))

    test_tsvpro_out = load_tsvpro_results(tsvout_path)
    expected_output = load_tsvpro_results((expected_tsvpro_outdir / f"{member_db}.tsvpro.tsv"))

    compare_tsvpro_output(test_tsvpro_out, )
    tsvout_path.unlink()
