import os
import json
import subprocess


def get_current_output(test_output_dir: str, current_output_path: str, input_path: str, applications: str, disable_precalc: bool) -> dict:
    disable_precalc = "--disable_precalc" if disable_precalc else ""
    command = f"nextflow run interproscan.nf --input {input_path} --applications {applications} {disable_precalc} " \
              f"--formats json --outdir {test_output_dir} --goterms --pathways -profile docker"
    if os.path.exists(str(current_output_path) + ".json"):
        os.remove(str(current_output_path) + ".json")
    subprocess.run(command, shell=True)
    with open(str(current_output_path) + ".json", 'r') as f:
        return json.load(f)


def get_expected_output(expected_output_path: str) -> dict:
    with open(str(expected_output_path) + ".json", 'r') as f:
        return json.load(f)


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


def compare(expected, current, ignore_elements: list, print_md5: bool, print_acc: bool, file=None):
    output = print if file is None else file.write
    for key in expected:
        if key == "md5" and print_md5:
            output(f"md5: {expected[key]}\n")
            print_md5 = False
        if key == "accession" and print_acc:
            output(f"accession: {expected[key]}\n")
            print_acc = False
        if key in ignore_elements:
            continue
        if key not in current:
            output(f"MISMATCH: Key '{key}' missing in current dict\n")
            print_md5 = True
            print_acc = True
            continue
        if isinstance(expected[key], dict):
            compare(expected[key], current[key], ignore_elements, print_md5, print_acc, file)
        elif isinstance(expected[key], list):
            if len(expected[key]) != len(current[key]):
                output((
                    f"MISMATCH: list length for key '{key}'\n"
                    f"expected: {len(expected[key])}\n"
                    f"current: {len(current[key])}\n"
                ))
            else:
                for i in range(len(expected[key])):
                    compare(expected[key][i], current[key][i], ignore_elements, print_md5, print_acc, file)
        else:
            if str(expected[key]).lower() != str(current[key]).lower():
                output(f"MISMATCH: for key '{key}'\n")
                output(f"  expected: {expected[key]}\n")
                output(f"  current: {current[key]}\n")
                print_md5 = True
                print_acc = True


def test_json_output(test_output_dir, input_path, expected_output_path, output_path, applications, disable_precalc):
    expected_output = get_expected_output(expected_output_path)
    current_output = get_current_output(test_output_dir, output_path, input_path, applications, disable_precalc)

    expected = json2dict(expected_output)
    current = json2dict(current_output)

    with open('tests/integration_tests/temp_expected.json', 'w') as file:
        json.dump(expected, file, indent=2)
    with open('tests/integration_tests/temp_current.json', 'w') as file:
        json.dump(current, file, indent=2)

    ignore_elements = []
    with open('tests/integration_tests/json_results.txt', 'w') as file:
        compare(expected, current, ignore_elements, False, False, file)

    assert expected == current
