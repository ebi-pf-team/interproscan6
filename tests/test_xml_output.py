import pytest
import subprocess
import xml.etree.ElementTree as ET
import os


def run_interproscan(input_file, output_file):
    command = f"nextflow run interproscan.nf --input {input_file} --applications antifam --formats xml --disable_precalc --output {output_file}"
    subprocess.run(command, shell=True)


@pytest.fixture
def output_data():
    output_file = "tests/test_outputs/xml_output"
    input_file = "tests/test_outputs/xml_test.fasta"
    if os.path.exists(output_file + ".xml"):
        os.remove(output_file + ".xml")
    run_interproscan(input_file, output_file)
    return parse_xml("tests/test_outputs/xml_output.xml")


@pytest.fixture
def expected_output_data():
    return parse_xml("tests/test_outputs/xml_output_expected.xml")


def parse_xml(xml_file):
    with open(xml_file, 'r') as f:
        tree = ET.parse(f)
    root = tree.getroot()
    data = {}

    ns = {"iprscan": "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6/schemas"}

    for protein in root.findall(".//iprscan:protein", ns):
        protein_id = protein.find(".//iprscan:xref", ns).attrib["id"]
        matches = []
        for match in protein.findall(".//iprscan:matches/iprscan:hmmer3-match", ns):
            match_data = {
                "sequence": protein.find(".//iprscan:sequence", ns).text,
                "xref_id": protein.find(".//iprscan:xref", ns).attrib["id"],
                "hmmer_match": []
            }
            match_data["hmmer_match"].append(match.find("iprscan:model-ac", ns).text)
            matches.append(match_data)
        data[protein_id] = matches
    return data


def compare_xml(output_data, expected_output_data):
    for protein_id, output_matches in output_data.items():
        assert protein_id in expected_output_data, f"Protein ID {protein_id} not found in expected output data"
        expected_matches = expected_output_data[protein_id]
        assert len(output_matches) == len(expected_matches), f"Protein {protein_id} has different number of matches"

        for output_match, expected_match in zip(output_matches, expected_matches):
            assert output_match.keys() == expected_match.keys(), f"Keys for protein {protein_id} do not match"
            for key in output_match.keys():
                if isinstance(output_match[key], list):
                    assert sorted(output_match[key]) == sorted(expected_match[key]), f"List elements for protein {protein_id} do not match"
                else:
                    assert output_match[key] == expected_match[key], f"Values for protein {protein_id} do not match"


def test_output(output_data, expected_output_data):
    compare_xml(output_data, expected_output_data)
