import argparse
import json
import csv
import os
import xml.etree.ElementTree as ET


def xml_output(matches: str, output_path: str):
    # A lot of changes! Need to be recreated!
    pass


def tsv_output(matches: str, output_path: str):
    tsv_output = output_path + '.tsv'

    with open(matches, 'r') as file:
        data = file.readlines()
        with open(tsv_output, 'w') as tsv_file:
            for member in data:
                member_matches = json.loads(member)
                for match in member_matches:
                    tsv_output = csv.writer(tsv_file, delimiter='\t')
                    output_line = match["tbl"] + match["sequence"]
                    # tsv_file.write(output_line)
                    tsv_output.writerow(output_line)


def json_output(matches: str, output_path: str):
    with open(matches, 'r') as file:
        data = file.read()

    output_with_format = output_path + '.json'
    # concatenated_data = {"interproscan-version": "6.0.0", 'results': data}

    with open(output_with_format, 'w') as json_file:
        json.dump(data, json_file, indent=2)


def gff3_output(matches: str, output_path: str):
    pass


def write_results(matches: str, output_format: str, output_path: str):
    output_format = output_format.upper()
    if "TSV" in output_format:
        tsv_output(matches, output_path)
    if "XML" in output_format:
        xml_output(matches, output_path)
    if "JSON" in output_format:
        json_output(matches, output_path)
    if "GFF3" in output_format:
        gff3_output(matches, output_path)


def main():
    parser = argparse.ArgumentParser(description="Write result file")
    parser.add_argument(
        "-results", "--results", type=str, help="matches result parsed"
    )
    parser.add_argument("-format", "--format", type=str, help="output format")
    parser.add_argument("-output_path", "--output_path", type=str, help="output path")

    args = parser.parse_args()
    write_results(args.results, args.format, args.output_path)


if __name__ == "__main__":
    main()
