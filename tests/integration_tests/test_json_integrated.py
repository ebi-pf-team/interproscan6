import json
import os
import shutil
import zipfile
import subprocess
from collections import OrderedDict
from pathlib import Path


def get_current_output(test_output_dir: str,
                       current_output_path: str,
                       input_path: str,
                       applications: str,
                       data_dir: str,
                       disable_precalc: bool) -> dict:
    disable_precalc = "--disable_precalc" if disable_precalc else ""
    command = f"nextflow run main.nf --input {input_path} --applications {applications} " \
              f"--formats json --outdir {test_output_dir} -profile docker \
              --datadir {data_dir} -resume"
    # if os.path.exists(str(current_output_path) + ".json"):
    #     os.remove(str(current_output_path) + ".json")
    # subprocess.run(command, shell=True)
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
    signature_fields = ["accession", "entry"]  # "name", "description"
    locations_fields = ["start", "end", "representative", "hmmStart", "hmmEnd", "hmmLength", "hmmBounds",
                        "evalue", "score", "envelopeStart", "envelopeEnd", "bias", "cigarAlignment", "level",
                        "sequenceFeature", "pvalue", "motifNumber", "queryAlignment", "targetAlignment"]

    result = {}
    for result_item in data.get("results", []):
        for match in result_item.get("matches", []):
            signature = match.get("signature", {})
            library = signature.get("signatureLibraryRelease", {}).get("library", "").lower().replace("-", "").replace(" ", "").replace("_", "")
            if not library:
                continue
            if library not in result:
                result[library] = {}

            accession = signature.get("accession")
            evalue = match.get("evalue", -1)
            score = match.get("score", -1)
            signature_filtered = OrderedDict((k, signature[k]) for k in signature_fields if k in signature)
            for xref in result_item.get("xref", []):
                xref_id = xref.get("id")
                if not xref_id:
                    continue
                if xref_id not in result[library]:
                    result[library][xref_id] = {}
                if accession not in result[library][xref_id]:
                    result[library][xref_id][accession] = {"locations": {}}
                for loc in match.get("locations", []):
                    start = loc.get("start", -1)
                    end = loc.get("end", -1)
                    location_key = f"s{start}-e{end}"
                    location_data = OrderedDict()
                    for field in locations_fields:
                        if field in loc:
                            location_data[field] = loc[field]
                    result[library][xref_id][accession]["locations"][location_key] = OrderedDict((k, loc[k]) for k in locations_fields if k in loc)
                    if "sites" in loc:
                        if "sites" not in result[library][xref_id][accession]:
                            result[library][xref_id][accession]["sites"] = {}
                        loc["sites"] = sorted(
                            loc["sites"],
                            key=lambda site: (
                                site.get("hmmStart", -1),
                                site.get("hmmEnd", -1)
                            )
                        )
                        result[library][xref_id][accession]["sites"][location_key] = loc["sites"]
                    if "location-fragments" in loc:
                        if "fragments" not in result[library][xref_id][accession]:
                            result[library][xref_id][accession]["fragments"] = {}
                        result[library][xref_id][accession]["fragments"][location_key] = loc["location-fragments"]
                result[library][xref_id][accession]["signature"] = signature_filtered
                result[library][xref_id][accession]["evalue"] = evalue
                result[library][xref_id][accession]["score"] = score
    return OrderedDict(sorted(result.items()))


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

                        # EVALUE (MATCH) CHECK
                        evalue1 = expected[library][xref_id][accession]["evalue"]
                        evalue2 = current[library][xref_id][accession]["evalue"]
                        if evalue1 != evalue2:
                            file.write(
                                f"MISMATCH ON MATCH EVALUE: {library} -> {xref_id} -> {accession}:\n\t Expected: {evalue1} \n\tCurrent : {evalue2}\n")
                        # SCORE (MATCH) CHECK
                        score1 = expected[library][xref_id][accession]["score"]
                        score2 = current[library][xref_id][accession]["score"]
                        if score1 != score2:
                            file.write(
                                f"MISMATCH ON MATCH SCORE: {library} -> {xref_id} -> {accession}:\n\tExpected: {score1} \n\tCurrent : {score2}\n")
                        # SIGNATURE CHECK
                        signature1 = expected[library][xref_id][accession]["signature"]
                        signature2 = current[library][xref_id][accession]["signature"]
                        if signature1 != signature2:
                            file.write(
                                f"MISMATCH ON SIGNATURE: {library} -> {xref_id} -> {accession}:\n\t Expected: {signature1} \n\t Current : {signature2}\n")
                        # LOCATIONS CHECK
                        locations1 = expected[library][xref_id][accession]["locations"]
                        locations2 = current[library][xref_id][accession]["locations"]
                        if locations1 != locations2:
                            file.write(f"MISMATCH ON LOCATIONS: {library} -> {xref_id} -> {accession}:\n\t Expected: {locations1} \n\t Current : {locations2}\n")

                        # SITES CHECK
                        sites1 = expected[library][xref_id][accession].get("sites", {})
                        sites2 = current[library][xref_id][accession].get("sites", {})
                        if sites1 != sites2:
                            file.write(f"MISMATCH ON SITES: {library} -> {xref_id} -> {accession}:\n\t Expected: {sites1} \n\t Current : {sites2}\n")

                        # FRAGMENTS CHECK
                        fragments1 = expected[library][xref_id][accession].get("fragments", {})
                        fragments2 = current[library][xref_id][accession].get("fragments", {})
                        if fragments1 != fragments2:
                            file.write(f"MISMATCH ON FRAGMENTS: {library} -> {xref_id} -> {accession}:\n\t Expected: {fragments1} \n\t Current : {fragments2}\n")

            # file.write(f"Skipping applications: {', '.join(skipping_libs)}")  # JUST TO DEBUG


def test_json_output(test_output_dir, input_path, expected_result_path, output_path, applications, data_dir, disable_precalc):
    expected_result = get_expected_result(expected_result_path)
    current_result = get_current_output(test_output_dir, output_path, input_path, applications, data_dir, disable_precalc)

    expected = json2dict(expected_result)
    current = json2dict(current_result)

    output_file = "tests/integration_tests/mismatches.txt"
    compare(expected, current, output_file)

    assert os.stat(output_file).st_size == 0
