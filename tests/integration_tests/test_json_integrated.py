import json
import subprocess


def get_current_output(current_output_path: str, input_path: str, applications: str, disable_precalc: bool) -> dict:
    disable_precalc = "--disable_precalc" if disable_precalc else ""
    command = f"nextflow run interproscan.nf --input {input_path} --applications {applications} {disable_precalc} --formats json,tsv-pro --output {current_output_path} --goterms --pathways"
    if os.path.exists(current_output_path + ".json"):
        os.remove(current_output_path + ".json")
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
        # Ordena todos os itens da lista recursivamente
        sorted_list = sorted([json2dict(item) for item in obj], key=lambda x: json.dumps(x))
        # Verifica se algum item tem a chave 'md5' e reordena com base nela se necessário
        if sorted_list and isinstance(sorted_list[0], dict) and 'md5' in sorted_list[0]:
            sorted_list = sorted(sorted_list, key=lambda x: x.get('md5', ''))
        return sorted_list
    else:
        return obj


def compare(expected, current, ignore_elements: list):
    for key in expected:
        if key in ignore_elements:
            continue
        if key not in current:
            print(f"Key '{key}' missing in current dict")
            continue
        if isinstance(expected[key], dict):
            compare(expected[key], current[key], ignore_elements)
        elif isinstance(expected[key], list):
            if len(expected[key]) != len(current[key]):
                print(f"List length mismatch for key '{key}'")
            else:
                for i in range(len(expected[key])):
                    compare(expected[key][i], current[key][i], ignore_elements)
        else:
            if str(expected[key]).lower() != str(current[key]).lower():
                print(f"Value mismatch for key '{key}'")
                print(f"  expected: {expected[key]}")
                print(f"  current: {current[key]}")


def test_json_output(input_path, expected_output_path, current_output_path, applications, disable_precalc):
    expected_output = get_expected_output(expected_output_path)
    current_output = get_current_output(current_output_path, input_path, applications, disable_precalc)

    expected = json2dict(expected_output)
    current = json2dict(current_output)

    ignore_elements = ['representative']
    print("Missing elements in current output:")
    compare(expected, current, ignore_elements)
    print("Extra elements in current output:")
    compare(current, expected, ignore_elements)

    with open('temp_expected.json', 'w') as file:
        json.dump(expected, file, indent=2)
    with open('temp_current.json', 'w') as file:
        json.dump(current, file, indent=2)

    assert expected == current  # Uncomment this line when output totally implemented
