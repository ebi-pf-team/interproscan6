import subprocess
import xml.etree.ElementTree as ET


def get_current_output(current_output_path: str, input_path: str, applications: str, disable_precalc: bool) -> ET.Element:
    disable_precalc = "--disable_precalc" if disable_precalc else ""
    command = f"nextflow run interproscan.nf --input {input_path} --applications {applications} {disable_precalc} --formats xml --output {current_output_path}"
    subprocess.run(command, shell=True)
    with open(str(current_output_path) + ".xml", 'r') as f:
        tree = ET.parse(f)
    return tree.getroot()


def get_expected_output(expected_output_path: str) -> ET.Element:
    with open(str(expected_output_path) + ".xml", 'r') as f:
        tree = ET.parse(f)
    return tree.getroot()


def xml2dict(element: ET.Element):
    result = {}
    if element.attrib:
        result.update({'attributes': element.attrib})
    for child in element:
        if child.tag not in result:
            result[child.tag] = []
        result[child.tag].append(xml2dict(child))
    for key, value in result.items():
        if isinstance(value, list) and all(isinstance(item, dict) for item in value):
            result[key] = sorted(value, key=lambda x: str(x))
    if 'locations' in result:
        result['locations'] = sorted(result['locations'], key=lambda x: int(x['hmmer3-location'][0]['attributes']['env-end']))
    return result


def compare(expected, current, ignore_elements: list):
    for key in expected:
        if key in ignore_elements:
            continue
        if key not in current:
            raise AssertionError(f"Key '{key}' missing in current dict")
        if isinstance(expected[key], dict):
            compare(expected[key], current[key], ignore_elements)
        elif isinstance(expected[key], list):
            if len(expected[key]) != len(current[key]):
                raise AssertionError(f"List length mismatch for key '{key}'")
            else:
                for i in range(len(expected[key])):
                    compare(expected[key][i], current[key][i], ignore_elements)
        else:
            if str(expected[key]).lower() != str(current[key]).lower():
                print(f"  expected: {expected[key]}")
                print(f"  current: {current[key]}")
                raise AssertionError(f"Value mismatch for key '{key}'")


def remove_namespace(element):
    for elem in element.iter():
        if '}' in elem.tag:
            elem.tag = elem.tag.split('}', 1)[1]
        elem.attrib = {k.split('}', 1)[-1].lower(): v for k, v in elem.attrib.items()}


def test_xml_output(input_path, expected_output_path, current_output_path, applications, disable_precalc):
    expected_output = get_expected_output(expected_output_path)
    current_output = get_current_output(current_output_path, input_path, applications, disable_precalc)
    remove_namespace(expected_output)
    remove_namespace(current_output)

    expected = xml2dict(expected_output)
    current = xml2dict(current_output)

    ignore_elements = ['representative', 'hmmer3-location-fragment', 'hmm-bounds', 'evalue']
    print("Missing elements in current output:")
    compare(expected, current, ignore_elements)
    print("Extra elements in current output:")
    compare(current, expected, ignore_elements)

    # assert expected == current  # Uncomment this line when output totally implemented
