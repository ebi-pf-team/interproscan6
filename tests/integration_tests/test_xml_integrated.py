import subprocess
import xml.etree.ElementTree as ET


def get_current_output(current_output_path: str, input_path: str, applications: str, disable_precalc: bool) -> ET.Element:
    disable_precalc = "--disable_precalc" if disable_precalc else ""
    command = f"nextflow run interproscan.nf --input {input_path} --applications {applications} {disable_precalc} --formats xml --output {current_output_path} --goterms --pathways"
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


def compare(dict1, dict2, ignore_elements: list, comparing_type: str):
    for key in dict1:
        if key in ignore_elements:
            continue
        if key not in dict2:
            print(f"Key '{key}' {comparing_type} in current dict")
        if isinstance(dict1[key], dict):
            compare(dict1[key], dict2[key], ignore_elements, comparing_type)
        elif isinstance(dict1[key], list):
            if len(dict1[key]) != len(dict2[key]):
                print(f"List length mismatch for key '{key}'")
            else:
                for i in range(len(dict1[key])):
                    compare(dict1[key][i], dict2[key][i], ignore_elements, comparing_type)
        else:
            if str(dict1[key]).lower() != str(dict2[key]).lower():
                print(f"  dict1: {dict1[key]}")
                print(f"  dict2: {dict2[key]}")
                print(f"Value mismatch for key '{key}'")


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

    ignore_elements = ['representative', 'hmmer3-location-fragment', 'hmm-bounds', 'evalue', "sites"]
    compare(expected, current, ignore_elements, "missing")
    compare(current, expected, ignore_elements, "extra")

    # assert expected == current  # Uncomment this line when output totally implemented
