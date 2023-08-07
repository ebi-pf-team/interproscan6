import csv
import argparse
import json
import xml.etree.ElementTree as ET

SERIAL_GROUP = "PROTEIN"
# Default for protein sequences are TSV, XML and GFF3, for nucleotide sequences GFF3 and XML.


def xml_output(members_info):
    # TODO filter xml infos
    with open(f"lookup.xml", 'w') as file:
        root = ET.Element('results')
        for result in members_info:
            elemento = ET.SubElement(root, 'result')
            elemento.text = result
        tree = ET.ElementTree(root)
        tree.write(file)


def tsv_output(members_info):
    # TODO filter tsv infos
    with open("output.tsv", 'w') as file:
        for entry in members_info:
            writer = csv.writer(file, delimiter='\t')
            writer.writerow(entry.values())
            writer.writerow("\n")


def json_output(members_info):
    # TODO filter json infos
    with open("output.json", 'w') as file:
        json.dump(members_info, file)


def gff3_output(members_info):
    pass


def write_results(members_info: list[dict], output_format: list):
    if len(output_format) == 0:
        if SERIAL_GROUP == "PROTEIN":
            output_format = ["TSV", "XML", "GFF3"]
        else:
            output_format = ["XML", "GFF3"]
    else:
        output_format = list(map(lambda x: x.upper(), output_format))

    if 'TSV' in output_format:
        tsv_output(members_info)
    if 'XML' in output_format:
        xml_output(members_info)
    if 'JSON' in output_format:
        json_output(members_info)
    if 'GFF3' in output_format:
        gff3_output(members_info)


def main():
    parser = argparse.ArgumentParser(description="Write result file(s)")
    parser.add_argument("results_txt", metavar="results_txt", type=str, help="result file")
    parser.add_argument(
        "output_formats", metavar="output_formats", type=str, help="the output formats"
    )

    args = parser.parse_args()
    write_results(args.results_txt, args.output_formats)


if __name__ == '__main__':
    main()
