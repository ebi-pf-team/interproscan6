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
            if seq_id not in results:
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
def expected_tsvpro_outdir(test_output_dir):
    return test_output_dir / "format_writer/tsv-pro"


@pytest.fixture
def tsvout_path(test_output_dir):
    dirpath = test_output_dir / "temp"
    filepath = dirpath / "tsv.unittest.ips6.tsv"
    dirpath.mkdir(parents=True, exist_ok=True)
    return filepath


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



def test_tsvpro_antifam_output(tsv_seq_match_dir, tsvout_path, expected_tsvpro_outdir):
    member_db = "ANTIFAM"
    seq_matches = load_seq_matches_dict(member_db, tsv_seq_match_dir)
    tsv_output.tsv_pro_output(seq_matches, tsvout_path)

    test_tsvpro_out = load_tsvpro_results(tsvout_path)
    expected_output = load_tsvpro_results((expected_tsvpro_outdir / f"{member_db}.tsvpro.tsv"))

    compare_tsvpro_output(test_tsvpro_out, expected_output)
    tsvout_path.unlink()


def test_tsvpro_cdd_output(tsv_seq_match_dir, tsvout_path, expected_tsvpro_outdir):
    member_db = "CDD"
    seq_matches = load_seq_matches_dict(member_db, tsv_seq_match_dir)
    tsv_output.tsv_pro_output(seq_matches, tsvout_path)

    test_tsvpro_out = load_tsvpro_results(tsvout_path)
    expected_output = load_tsvpro_results((expected_tsvpro_outdir / f"{member_db}.tsvpro.tsv"))

    compare_tsvpro_output(test_tsvpro_out, expected_output)
    tsvout_path.unlink()


def test_tsvpro_coils_output(tsv_seq_match_dir, tsvout_path, expected_tsvpro_outdir):
    member_db = "COILS"
    seq_matches = load_seq_matches_dict(member_db, tsv_seq_match_dir)
    tsv_output.tsv_pro_output(seq_matches, tsvout_path)

    test_tsvpro_out = load_tsvpro_results(tsvout_path)
    expected_output = load_tsvpro_results((expected_tsvpro_outdir / f"{member_db}.tsvpro.tsv"))

    compare_tsvpro_output(test_tsvpro_out, expected_output)
    tsvout_path.unlink()


def test_tsvpro_funfam_output(tsv_seq_match_dir, tsvout_path, expected_tsvpro_outdir):
    member_db = "FUNFAM"
    seq_matches = load_seq_matches_dict(member_db, tsv_seq_match_dir)
    tsv_output.tsv_pro_output(seq_matches, tsvout_path)

    test_tsvpro_out = load_tsvpro_results(tsvout_path)
    expected_output = load_tsvpro_results((expected_tsvpro_outdir / f"{member_db}.tsvpro.tsv"))

    compare_tsvpro_output(test_tsvpro_out, expected_output)
    tsvout_path.unlink()


def test_tsvpro_gene3d_output(tsv_seq_match_dir, tsvout_path, expected_tsvpro_outdir):
    member_db = "GENE3D"
    seq_matches = load_seq_matches_dict(member_db, tsv_seq_match_dir)
    tsv_output.tsv_pro_output(seq_matches, tsvout_path)

    test_tsvpro_out = load_tsvpro_results(tsvout_path)
    expected_output = load_tsvpro_results((expected_tsvpro_outdir / f"{member_db}.tsvpro.tsv"))

    compare_tsvpro_output(test_tsvpro_out, expected_output)
    tsvout_path.unlink()


def test_tsvpro_hamap_output(tsv_seq_match_dir, tsvout_path, expected_tsvpro_outdir):
    member_db = "HAMAP"
    seq_matches = load_seq_matches_dict(member_db, tsv_seq_match_dir)
    tsv_output.tsv_pro_output(seq_matches, tsvout_path)

    test_tsvpro_out = load_tsvpro_results(tsvout_path)
    expected_output = load_tsvpro_results((expected_tsvpro_outdir / f"{member_db}.tsvpro.tsv"))

    compare_tsvpro_output(test_tsvpro_out, expected_output)
    tsvout_path.unlink()


def test_tsvpro_mobidb_output(tsv_seq_match_dir, tsvout_path, expected_tsvpro_outdir):
    member_db = "MOBIDB_LITE"
    seq_matches = load_seq_matches_dict(member_db, tsv_seq_match_dir)
    tsv_output.tsv_pro_output(seq_matches, tsvout_path)

    test_tsvpro_out = load_tsvpro_results(tsvout_path)
    expected_output = load_tsvpro_results((expected_tsvpro_outdir / f"{member_db}.tsvpro.tsv"))

    compare_tsvpro_output(test_tsvpro_out, expected_output)
    tsvout_path.unlink()


def test_tsvpro_ncbifam_output(tsv_seq_match_dir, tsvout_path, expected_tsvpro_outdir):
    member_db = "NCBIFAM"
    seq_matches = load_seq_matches_dict(member_db, tsv_seq_match_dir)
    tsv_output.tsv_pro_output(seq_matches, tsvout_path)

    test_tsvpro_out = load_tsvpro_results(tsvout_path)
    expected_output = load_tsvpro_results((expected_tsvpro_outdir / f"{member_db}.tsvpro.tsv"))

    compare_tsvpro_output(test_tsvpro_out, expected_output)
    tsvout_path.unlink()


def test_tsvpro_panther_output(tsv_seq_match_dir, tsvout_path, expected_tsvpro_outdir):
    member_db = "PANTHER"
    seq_matches = load_seq_matches_dict(member_db, tsv_seq_match_dir)
    tsv_output.tsv_pro_output(seq_matches, tsvout_path)

    test_tsvpro_out = load_tsvpro_results(tsvout_path)
    expected_output = load_tsvpro_results((expected_tsvpro_outdir / f"{member_db}.tsvpro.tsv"))

    compare_tsvpro_output(test_tsvpro_out, expected_output)
    tsvout_path.unlink()


def test_tsvpro_pfam_output(tsv_seq_match_dir, tsvout_path, expected_tsvpro_outdir):
    member_db = "PFAM"
    seq_matches = load_seq_matches_dict(member_db, tsv_seq_match_dir)
    tsv_output.tsv_pro_output(seq_matches, tsvout_path)

    test_tsvpro_out = load_tsvpro_results(tsvout_path)
    expected_output = load_tsvpro_results((expected_tsvpro_outdir / f"{member_db}.tsvpro.tsv"))

    compare_tsvpro_output(test_tsvpro_out, expected_output)
    tsvout_path.unlink()


def test_tsvpro_phobius_output(tsv_seq_match_dir, tsvout_path, expected_tsvpro_outdir):
    member_db = "PHOBIUS"
    seq_matches = load_seq_matches_dict(member_db, tsv_seq_match_dir)
    tsv_output.tsv_pro_output(seq_matches, tsvout_path)

    test_tsvpro_out = load_tsvpro_results(tsvout_path)
    expected_output = load_tsvpro_results((expected_tsvpro_outdir / f"{member_db}.tsvpro.tsv"))

    compare_tsvpro_output(test_tsvpro_out, expected_output)
    tsvout_path.unlink()


def test_tsvpro_pirsf_output(tsv_seq_match_dir, tsvout_path, expected_tsvpro_outdir):
    member_db = "PIRSF"
    seq_matches = load_seq_matches_dict(member_db, tsv_seq_match_dir)
    tsv_output.tsv_pro_output(seq_matches, tsvout_path)

    test_tsvpro_out = load_tsvpro_results(tsvout_path)
    expected_output = load_tsvpro_results((expected_tsvpro_outdir / f"{member_db}.tsvpro.tsv"))

    compare_tsvpro_output(test_tsvpro_out, expected_output)
    tsvout_path.unlink()


def test_tsvpro_pirsr_output(tsv_seq_match_dir, tsvout_path, expected_tsvpro_outdir):
    member_db = "PIRSR"
    seq_matches = load_seq_matches_dict(member_db, tsv_seq_match_dir)
    tsv_output.tsv_pro_output(seq_matches, tsvout_path)

    test_tsvpro_out = load_tsvpro_results(tsvout_path)
    expected_output = load_tsvpro_results((expected_tsvpro_outdir / f"{member_db}.tsvpro.tsv"))

    compare_tsvpro_output(test_tsvpro_out, expected_output)
    tsvout_path.unlink()


def test_tsvpro_prints_output(tsv_seq_match_dir, tsvout_path, expected_tsvpro_outdir):
    member_db = "PRINTS"
    seq_matches = load_seq_matches_dict(member_db, tsv_seq_match_dir)
    tsv_output.tsv_pro_output(seq_matches, tsvout_path)

    test_tsvpro_out = load_tsvpro_results(tsvout_path)
    expected_output = load_tsvpro_results((expected_tsvpro_outdir / f"{member_db}.tsvpro.tsv"))

    compare_tsvpro_output(test_tsvpro_out, expected_output)
    tsvout_path.unlink()


def test_tsvpro_prosite_patterns_output(tsv_seq_match_dir, tsvout_path, expected_tsvpro_outdir):
    member_db = "PROSITE_PATTERNS"
    seq_matches = load_seq_matches_dict(member_db, tsv_seq_match_dir)
    tsv_output.tsv_pro_output(seq_matches, tsvout_path)

    test_tsvpro_out = load_tsvpro_results(tsvout_path)
    expected_output = load_tsvpro_results((expected_tsvpro_outdir / f"{member_db}.tsvpro.tsv"))

    compare_tsvpro_output(test_tsvpro_out, expected_output)
    tsvout_path.unlink()


def test_tsvpro_prosite_profiles_output(tsv_seq_match_dir, tsvout_path, expected_tsvpro_outdir):
    member_db = "PROSITE_PROFILES"
    seq_matches = load_seq_matches_dict(member_db, tsv_seq_match_dir)
    tsv_output.tsv_pro_output(seq_matches, tsvout_path)

    test_tsvpro_out = load_tsvpro_results(tsvout_path)
    expected_output = load_tsvpro_results((expected_tsvpro_outdir / f"{member_db}.tsvpro.tsv"))

    compare_tsvpro_output(test_tsvpro_out, expected_output)
    tsvout_path.unlink()


def test_tsvpro_sfld_output(tsv_seq_match_dir, tsvout_path, expected_tsvpro_outdir):
    member_db = "SFLD"
    seq_matches = load_seq_matches_dict(member_db, tsv_seq_match_dir)
    tsv_output.tsv_pro_output(seq_matches, tsvout_path)

    test_tsvpro_out = load_tsvpro_results(tsvout_path)
    expected_output = load_tsvpro_results((expected_tsvpro_outdir / f"{member_db}.tsvpro.tsv"))

    compare_tsvpro_output(test_tsvpro_out, expected_output)
    tsvout_path.unlink()


def test_tsvpro_signalp_output(tsv_seq_match_dir, tsvout_path, expected_tsvpro_outdir):
    member_db = "SIGNALP"
    seq_matches = load_seq_matches_dict(member_db, tsv_seq_match_dir)
    tsv_output.tsv_pro_output(seq_matches, tsvout_path)

    test_tsvpro_out = load_tsvpro_results(tsvout_path)
    expected_output = load_tsvpro_results((expected_tsvpro_outdir / f"{member_db}.tsvpro.tsv"))

    compare_tsvpro_output(test_tsvpro_out, expected_output)
    tsvout_path.unlink()


def test_tsvpro_smart_output(tsv_seq_match_dir, tsvout_path, expected_tsvpro_outdir):
    member_db = "SMART"
    seq_matches = load_seq_matches_dict(member_db, tsv_seq_match_dir)
    tsv_output.tsv_pro_output(seq_matches, tsvout_path)

    test_tsvpro_out = load_tsvpro_results(tsvout_path)
    expected_output = load_tsvpro_results((expected_tsvpro_outdir / f"{member_db}.tsvpro.tsv"))

    compare_tsvpro_output(test_tsvpro_out, expected_output)
    tsvout_path.unlink()


def test_tsvpro_superfamily_output(tsv_seq_match_dir, tsvout_path, expected_tsvpro_outdir):
    member_db = "SUPERFAMILY"
    seq_matches = load_seq_matches_dict(member_db, tsv_seq_match_dir)
    tsv_output.tsv_pro_output(seq_matches, tsvout_path)

    test_tsvpro_out = load_tsvpro_results(tsvout_path)
    expected_output = load_tsvpro_results((expected_tsvpro_outdir / f"{member_db}.tsvpro.tsv"))

    compare_tsvpro_output(test_tsvpro_out, expected_output)
    tsvout_path.unlink()
