"""Test the python script that coordinates writing the output.

These test are intened to be run from the root of the repository using:
pytest -v
"""

import json

from pathlib import Path

import pytest

from interproscan.scripts.output.format_writer import tsv_output


def load_tsv_results(tsv_path: Path) -> dict[str, dict[str, str]]:
    results = {}
    with open(tsv_path, "r") as fh:
        for line in fh:
            data = line.split()
            seq_id = data[0]
            if seq_id not in results:
                # {memberdb: {sigacc: {'start-end': {match}}}}
                results[seq_id] = {}
            mem_db = data[3]
            if mem_db not in results[seq_id]:
                results[seq_id][mem_db] = {}
            sig_acc = data[4]
            if mem_db not in results[seq_id][mem_db]:
                results[seq_id][mem_db][sig_acc] = {}
            start, end = data[6], data[7]
            position = f"{start}-{end}"
            if position not in results[seq_id][mem_db][sig_acc]:
                results[seq_id][mem_db][sig_acc][position] = {
                    'start': start,
                    'end': end,
                    'md5': data[1],
                    'len': data[2],
                    'desc': data[4],
                    'score': data[8]
                }
    return results


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
            if seq_id not in results:
                # {memberdb: {sigacc: {'start-end': {match}}}}
                results[seq_id] = {}
            if mem_db not in results[seq_id]:
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
def tsvout_input(test_input_dir):
    _path = test_input_dir / "format_writer/protein/seq_match_dict.json"
    with open(_path, "r") as fh:
        _input = json.load(fh)
    return _input


@pytest.fixture
def tsvout_path(test_output_dir):
    dirpath = test_output_dir / "temp"
    filepath = dirpath / "tsv.unittest.ips6.tsv"
    dirpath.mkdir(parents=True, exist_ok=True)
    return filepath


@pytest.fixture
def expected_tsv_out(test_output_dir):
    _path = test_output_dir / "format_writer/expected_tsv.tsv"
    return load_tsv_results(_path)


@pytest.fixture
def expected_tsvpro_out(test_output_dir):
    _path = test_output_dir / "format_writer/expected_tsvpro.tsv"
    return load_tsvpro_results(_path)


def test_tsv_output(tsvout_input, tsvout_path, expected_tsv_out):
    tsv_output.tsv_output(tsvout_input, tsvout_path)

    test_tsv_output = load_tsv_results(tsvout_path)

    for seq_id, member_dbs in test_tsv_output.items():
        assert seq_id in expected_tsv_out, f"seq_id {seq_id} not in expected_tsv_out"
        for member_db, sig_accs in member_dbs.items():
            assert member_db in expected_tsv_out[seq_id], f"member_db {member_db} not in expected_tsv_out[{seq_id}]"
            for sig_acc, positions in sig_accs.items():
                assert sig_acc in expected_tsv_out[seq_id][member_db], f"sig_acc {sig_acc} not in expected_tsv_out[{seq_id}][{member_db}]"
                for position, match in positions.items():
                    expected_match = expected_tsv_out[seq_id][member_db][sig_acc][position]
                    assert position in expected_tsv_out[seq_id][member_db][sig_acc], f"position {position} not in expected_tsv_out[{seq_id}][{member_db}][{sig_acc}]"
                    assert all(
                        match[key] == expected_match[key]
                        for key in ['start', 'end', 'md5', 'len', 'desc', 'score']
                    ), f"Mismatch in match details for {seq_id}, {member_db}, {sig_acc}, {position}"

    tsvout_path.unlink()


def test_tsvpro_output(tsvout_input, tsvout_path, expected_tsvpro_out):
    tsv_output.tsv_pro_output(tsvout_input, tsvout_path)

    test_tsv_output = load_tsvpro_results(tsvout_path)

    for seq_id, member_dbs in test_tsv_output.items():
        assert seq_id in expected_tsvpro_out, f"seq_id {seq_id} not in expected_tsvpro_out"
        for member_db, sig_accs in member_dbs.items():
            assert member_db in expected_tsvpro_out[seq_id], f"member_db {member_db} not in expected_tsvpro_out[{seq_id}]"
            for sig_acc, positions in sig_accs.items():
                assert sig_acc in expected_tsvpro_out[seq_id][member_db], f"sig_acc {sig_acc} not in expected_tsvpro_out[{seq_id}][{member_db}]"
                for position, match in positions.items():
                    expected_match = expected_tsvpro_out[seq_id][member_db][sig_acc][position]
                    assert position in expected_tsvpro_out[seq_id][member_db][sig_acc], f"position {position} not in expected_tsvpro_out[{seq_id}][{member_db}][{sig_acc}]"
                    assert all(
                        match[key] == expected_match[key]
                        for key in ['release', 'fragments', 'score', 'hmm_start', 'hmm_end', 'hmm_len', 'location_score', 'location_evalue']
                    ), f"Mismatch in match details for {seq_id}, {member_db}, {sig_acc}, {position}"

    tsvout_path.unlink()
