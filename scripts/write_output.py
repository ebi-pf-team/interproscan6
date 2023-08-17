import argparse
import csv
import json
import xml.etree.ElementTree as ET

SERIAL_GROUP = "PROTEIN"
# Default for protein sequences are TSV, XML and GFF3, for nucleotide sequences GFF3 and XML.


def xml_output(matches_parsed):
    pass
    # it will change! will receive a dict with all infos and parse to each output format in this script


def tsv_output(matches_parsed):
    pass
    # it will change! will receive a dict with all infos and parse to each output format in this script


def json_output(members_info):
    pass  # focusing in xml and tsv for now (similar output info)


def gff3_output(members_info):
    pass  # focusing in xml and tsv for now (similar output info)


def write_results(members_info: list[dict], output_format: list):
    if len(output_format) > 0:
        output_format = list(map(lambda x: x.upper(), output_format))
    else:
        if SERIAL_GROUP == "PROTEIN":
            output_format = ["TSV", "XML", "GFF3"]
        else:
            output_format = ["XML", "GFF3"]

    if "TSV" in output_format:
        tsv_output(members_info)
    if "XML" in output_format:
        xml_output(members_info)
    if "JSON" in output_format:
        json_output(members_info)
    if "GFF3" in output_format:
        gff3_output(members_info)


def main():
    parser = argparse.ArgumentParser(description="Write result file(s)")
    parser.add_argument(
        "results", metavar="results", type=str, help="matches result parsed"
    )
    parser.add_argument(
        "formats", metavar="formats", type=str, help="output format(s)"
    )

    args = parser.parse_args()
    write_results(args.results_txt, args.output_formats)


if __name__ == "__main__":
    main()
