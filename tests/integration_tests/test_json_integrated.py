import json
import os
import shutil
import zipfile
import subprocess
from pathlib import Path
from collections import defaultdict


def get_current_output(test_output_dir: str,
                       current_output_path: str,
                       input_path: str,
                       applications: str,
                       data_dir: str,
                       disable_precalc: bool) -> dict:
    disable_precalc = "--disable_precalc" if disable_precalc else ""
    command = f"nextflow run main.nf --input {input_path} --applications {applications} " \
              f"--formats json --outdir {test_output_dir} -profile docker \
              --datadir {data_dir}"
    if os.path.exists(str(current_output_path) + ".json"):
        os.remove(str(current_output_path) + ".json")
    subprocess.run(command, shell=True)
    with open(str(current_output_path) + ".json", 'r') as f:
        return json.load(f)


def get_expected_result(expected_result_path: str) -> dict:
    path = Path(expected_result_path)
    if path.suffix == ".zip":
        with zipfile.ZipFile(path, 'r') as zip_ref:
            extracted_dir = path.parent / "extracted_temp"
            zip_ref.extractall(extracted_dir)
            json_files = [f for f in os.listdir(extracted_dir) if f.endswith(".json")]
            json_file = json_files[0]
        expected_result_file = os.path.join(extracted_dir, json_file)
    else:
        expected_result_file = expected_result_path + ".json"
    with open(expected_result_file, 'r') as f:
        expected_result = json.load(f)
    shutil.rmtree(extracted_dir, ignore_errors=True)

    return expected_result


def json2dict(data):
    result = {}
    for result_item in data.get("results", []):
        for match in result_item.get("matches", []):
            signature = match.get("signature", {})
            library = signature.get("signatureLibraryRelease", {}).get("library").lower().replace("-", "").replace(" ", "")
            name = signature.get("name", "")
            description = signature.get("description", "")
            entry = signature.get("entry", {})
            accession = signature.get("accession")
            if library is None or accession is None:
                continue
            if library not in result:
                result[library] = {}
            for xref in result_item.get("xref", []):
                xref_id = xref.get("id")
                if xref_id is None:
                    continue
                if xref_id not in result[library]:
                    result[library][xref_id] = {}
                locations = sorted(
                    match.get("locations", []),
                    key=lambda loc: (loc.get("start", float('inf')), loc.get("end", float('inf')))
                )
                result[library][xref_id][accession]["locations"] = locations
                result[library][xref_id][accession]["signature"] = signature
                result[library][xref_id][accession]["entry"] = entry
                result[library][xref_id][accession]["name"] = name
                result[library][xref_id][accession]["description"] = description

    return result


def compare(expected: dict,
            current: dict,
            output_file):

        skipping_libs = []
        all_libraries = set(expected.keys()).union(set(current.keys()))
        with open(output_file, "w") as file:
            for library in all_libraries:
                if library not in current:
                    skipping_libs.append(library)
                    continue
                all_xref_ids = set(expected[library].keys()).union(set(current[library].keys()))
                for xref_id in all_xref_ids:
                    if xref_id not in expected[library]:
                        file.write(f"EXTRA SEQUENCE: {library}: xref.id '{xref_id}'\n")
                        continue
                    if xref_id not in current[library]:
                        file.write(f"MISSING SEQUENCE: {library}: xref.id '{xref_id}'\n")
                        continue
                    all_accessions = set(expected[library][xref_id].keys()).union(set(current[library][xref_id].keys()))
                    for accession in all_accessions:
                        if accession not in expected[library][xref_id]:
                            file.write(f"EXTRA ACCESSION: {library} -> {xref_id}: accession '{accession}'\n")
                            continue
                        if accession not in current[library][xref_id]:
                            file.write(f"MISSING ACCESSION: {library} -> {xref_id}: accession '{accession}'\n")
                            continue
                        locations1 = expected[library][xref_id][accession].sort()
                        locations2 = current[library][xref_id][accession].sort()
                        if locations1 != locations2:
                            file.write(f"MISMATCH ON LOCATIONS: {library} -> {xref_id} -> {accession}:\n Expected: {locations1} \nCurrent.: {locations2}\n")
            file.write(f"Skipping applications: {', '.join(skipping_libs)}")


def test_json_output(test_output_dir, input_path, expected_result_path, output_path, applications, data_dir, disable_precalc):
    expected_result = get_expected_result(expected_result_path)
    current_result = get_current_output(test_output_dir, output_path, input_path, applications, data_dir, disable_precalc)

    expected = json2dict(expected_result)
    current = json2dict(current_result)

    output_file = "tests/integration_tests/mismatches.txt"
    compare(expected, current, output_file)

    assert os.stat(output_file).st_size == 0
