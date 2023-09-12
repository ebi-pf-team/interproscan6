import argparse
import json
import os
import csv
import xml.etree.ElementTree as ET


def xml_output(matches: list, output_path: str):
    # A lot of changes! Need to be recreated!
    pass

def tsv_output(matches: list, output_path: str):
    pass
    # A lot of changes! Need to be recreated!


def json_output(matches: list, output_path: str):
    output_with_format = output_path + '.json'
    with open(output_with_format, 'wb') as o:
        for match in matches:
            with open(match, "r") as m:
                o.write(json.load(m))


def gff3_output(matches: list, output_path: str):
    pass


def write_results(matches: list, output_format: str, output_path: str):
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
        "-results", "--results", nargs="*", help="matches result parsed"
    )
    parser.add_argument("-format", "--format", type=str, help="output format")
    parser.add_argument("-output_path", "--output_path", type=str, help="output path")

    args = parser.parse_args()
    write_results(args.results, args.format, args.output_path)


if __name__ == "__main__":
    main()
