import json
import os
import shutil
import zipfile
import subprocess
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


def json2dict(obj):
    if isinstance(obj, dict):
        result = {}
        for key, value in obj.items():
            result[key] = json2dict(value)
        return dict(sorted(result.items()))
    elif isinstance(obj, list):
        sorted_list = sorted([json2dict(item) for item in obj], key=lambda x: json.dumps(x))
        if sorted_list and isinstance(sorted_list[0], dict) and 'md5' in sorted_list[0]:
            sorted_list = sorted(sorted_list, key=lambda x: x.get('md5', ''))
        if sorted_list and isinstance(sorted_list[0], dict) and 'model-ac' in sorted_list[0]:
            sorted_list = sorted(sorted_list, key=lambda x: x.get('model-ac', ''))
        return sorted_list
    else:
        return obj


def compare(expected: dict,
            current: dict,
            ignore_elements: list,
            applications: list,
            seq_info=None,
            output_file=None):

    seq_info = seq_info or {}

    def update_seq_info(obj, seq_info):
        if "xref" in obj and isinstance(obj["xref"], list) and obj["xref"]:
            seq_info["xref_id"] = obj["xref"][0].get("id", "-")
        if "matches" in obj and isinstance(obj["matches"], list) and obj["matches"]:
            match = obj["matches"][0]
            seq_info["model_ac"] = match.get("model-ac", "-")
            if "signature" in match:
                seq_info["library"] = match["signature"].get("signatureLibraryRelease", {}).get("library", "-")

    def filter_matches(matches, applications):
        return [
            match for match in matches
            if match.get("signature", {}).get("signatureLibraryRelease", {}).get("library", "-").lower() in applications]

    update_seq_info(expected, seq_info)
    update_seq_info(current, seq_info)

    # just creating temp files for debugging
    with open('tests/integration_tests/temp_current.json', 'w') as file:
        json.dump(current, file, indent=2)
    with open('tests/integration_tests/temp_expected.json', 'w') as file:
        for key in expected:
            if key in ignore_elements:
                continue
            if key not in current:
                output_file.write(
                    f"{key}\tKey missing\t-\t{seq_info.get('xref_id', '-')}\t"
                    f"{seq_info.get('model_ac', '-')}\t{seq_info.get('library', '-')}\n"
                )
                continue

            if isinstance(expected[key], dict):
                compare(expected[key], current[key], ignore_elements, applications, seq_info, output_file)
            elif isinstance(expected[key], list):
                if key == "matches":
                    filtered_expected = filter_matches(expected[key], applications)
                    json.dump(filtered_expected, file, indent=2)

                    if len(filtered_expected) != len(current[key]):
                        output_file.write(
                            f"{key}\tList length mismatch\t-\t{seq_info.get('xref_id', '-')}\t"
                            f"{seq_info.get('model_ac', '-')}\t{seq_info.get('library', '-')}\n"
                        )
                    else:
                        for i in range(len(filtered_expected)):
                            compare(filtered_expected[i], current[key][i], ignore_elements, applications, seq_info, output_file)
                else:
                    if len(expected[key]) != len(current[key]):
                        output_file.write(
                            f"{key}\tList length mismatch\t-\t{seq_info.get('xref_id', '-')}\t"
                            f"{seq_info.get('model_ac', '-')}\t{seq_info.get('library', '-')}\n"
                        )
                    else:
                        for i in range(len(expected[key])):
                            compare(expected[key][i], current[key][i], ignore_elements, applications, seq_info, output_file)
            else:
                if str(expected[key]).lower().strip() != str(current[key]).lower().strip():
                    output_file.write(
                        f"{key}\t{expected[key]}\t{current[key]}\t{seq_info.get('xref_id', '-')}\t"
                        f"{seq_info.get('model_ac', '-')}\t{seq_info.get('library', '-')}\n"
                    )


def test_json_output(test_output_dir, input_path, expected_result_path, output_path, applications, data_dir, disable_precalc):
    expected_result = get_expected_result(expected_result_path)
    current_result = get_current_output(test_output_dir, output_path, input_path, applications, data_dir, disable_precalc)

    expected = json2dict(expected_result)
    current = json2dict(current_result)

    ignore_fields = ["postProcessed", 'dc-status', 'representative', 'library']
    applications_list = applications.split(',')
    with open("tests/integration_tests/mismatches.txt", "w") as file:
        compare(expected, current, ignore_fields, applications_list, output_file=file)

    assert expected == current
