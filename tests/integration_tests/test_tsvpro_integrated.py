"""Test the the final TSV-PRO output from IPS6 matches the expected output.

When testing against i5, all member databases were enabled except:
* TMHMM because IPS6 migrated from TMHMM to DeepTMHMM
* SignalP as IPS6 migrated from version 4 to version 6

These test are intended to be run from the root of the repository using:
python -m pytest -v
"""


import subprocess

from pathlib import Path

import pytest


def run_nextflow(
    outdir: str,
    input_path: str,
    applications: str,
    disable_precalc: bool
) -> dict:
    disable_precalc = "--disable_precalc" if disable_precalc else ""
    command = f"nextflow run interproscan.nf " \
              f"--input {input_path} --applications {applications} {disable_precalc} " \
              f"--formats tsv-pro " \
              "--goterms --pathways " \
              "-profile docker " \
              f"--outdir {str(outdir)}"
    with open("nextflow-run", "w") as fh:
        subprocess.run(command, shell=True, stdout=fh)


def load_tsvpro_results(tsv_path: Path) -> dict[str, dict[str, str]]:
    results = {}
    with open(tsv_path, "r") as fh:
        for line in fh:
            data = line.split()
            seq_id = data[3]
            mem_db = data[0]
            if mem_db == "PROSITEPATTERNS":
                mem_db = "PROSITE_PATTERNS"
            if mem_db == "PROSITEPROFILES":
                mem_db = "PROSITE_PROFILES"
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


def compare(
    expected:  dict[str, [dict[str, dict]]],
    current:  dict[str, [dict[str, dict]]],
    elements=['release', 'fragments', 'hmm_start', 'hmm_end', 'hmm_len']
):
    mismatch = False
    for seq_id, member_dbs in current.items():
        if seq_id not in expected:
            print(f"Seq id {seq_id} from the current output is not in the expected output")
            mismatch = True
            continue

        for member_db, sig_accs in member_dbs.items():
            if member_db not in expected[seq_id]:
                print(f"Member db {member_db} not in the expected output for {seq_id}")
                mismatch = True
                continue

            for sig_acc, positions in sig_accs.items():
                if sig_acc not in expected[seq_id][member_db]:
                    print(
                        f"Signature accession {sig_acc} not in the "
                        f"expected output for {seq_id}, {member_db}"
                    )
                    mismatch = True
                    continue

                for position, match in positions.items():
                    if position in expected[seq_id][member_db][sig_acc]:
                        expected_match = expected[seq_id][member_db][sig_acc][position]
                        if any(
                            match[key] != expected_match[key]
                            for key in elements
                        ):
                            print(
                                f"Mismatch in match details for {seq_id}, "
                                f"{member_db}, {sig_acc}, {position}"
                            )
                            mismatch = True
                        if match['score'] != expected_match['score']:
                            print(
                                f"Score differs between current and expected output "
                                f"for {seq_id}, {member_db}, {sig_acc} (position: {position}): "
                                f"expected: {expected_match['score']}, current: {match['score']}"
                            )
                        if match['location_score'] != expected_match['location_score']:
                            print(
                                f"Location score differs between current and expected output "
                                f"for {seq_id}, {member_db}, {sig_acc} (position: {position}): "
                                f"expected: {expected_match['score']}, current: {match['score']}"
                            )
                        if match['location_evalue'] != expected_match['location_evalue']:
                            print(
                                f"Location E-value differs between current and expected output "
                                f"for {seq_id}, {member_db}, {sig_acc} (position: {position}): "
                                f"expected: {expected_match['score']}, current: {match['score']}"
                            )
                    else:
                        print(
                            f"Position {position} not in expected output "
                            f"for {seq_id}, {member_db}, {sig_acc}"
                        )
                        mismatch = True

    for seq_id, member_dbs in expected.items():
        if seq_id not in current:
            print(
                f"Seq id {seq_id} is in the expected output "
                f"but is not in current output"
            )
            mismatch = True
            continue
        for member_db, sig_accs in member_dbs.items():
            if member_db not in current[seq_id]:
                print(
                    f"Member db {member_db} is in the expected "
                    f"output for {seq_id} but is not in the "
                    "current output"
                )
                mismatch = True
                continue
            for sig_acc, positions in sig_accs.items():
                if sig_acc not in current[seq_id][member_db]:
                    print(
                        f"Signature accession {sig_acc} is in the "
                        f"expected output for {seq_id}, {member_db} "
                        "but is not in the current output"
                    )
                    mismatch = True
                    continue
                for position, match in positions.items():
                    if position not in current[seq_id][member_db][sig_acc]:
                        print(
                            f"Position {position} is found in the expected output "
                            " but is not in the current output "
                            f"for {seq_id}, {member_db}, {sig_acc}"
                        )
                        mismatch = True

    return mismatch


def test_tsv_output(
    input_path,
    expected_output_path,
    current_output_path,
    applications,
    disable_precalc
):
    """Input parameters are defined in conf.py"""
    expected_output = load_tsvpro_results(Path(f"{expected_output_path}.tsvpro.tsv"))

    outdir = current_output_path.parent
    run_nextflow(
        outdir,
        input_path,
        applications,
        disable_precalc
    )
    output_path = outdir / f"{input_path.name}.ips6.tsv-pro.tsv"

    current_output = load_tsvpro_results(output_path)

    mismatch = compare(expected_output, current_output)
    if mismatch:
        pytest.fail("Current result did not match the expected result")
