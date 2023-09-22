import argparse
import json
import csv
import os
from datetime import datetime
import xml.etree.ElementTree as ET


def tsv_output(matches: str, output_path: str):
    tsv_output = os.path.join(output_path + '.tsv')
    with open(tsv_output, 'w') as tsv_file:
        current_date = datetime.now().strftime('%d-%m-%Y')
        for match in matches:
            seq_id = match["sequence"][0].split()[0]
            md5 = match["sequence"][2]
            seq_len = match["sequence"][3]
            for acc_match in match["acc_matches"]:
                acc = acc_match["accession"].split(".")[0]
                try:
                    signature_desc = acc_match["signature_desc"]
                    interpro_acc = acc_match["interpro_annotations_acc"]
                except:
                    signature_desc = "-"
                    interpro_acc = "-"
                for domain in acc_match["domains"]:
                    ali_from = domain["ali_from"]
                    ali_to = domain["ali_to"]
                    i_evalue = domain["iEvalue"]
                tsv_file.write(f"{seq_id}\t{md5}\t{seq_len}\t{acc}\t{signature_desc}\t{ali_from}\t{ali_to}\t{i_evalue}\t{current_date}\t{interpro_acc}\n")
        # tsv_output = csv.writer(tsv_file, delimiter='\t')
        # tsv_output.writerow(info)


def json_output(matches: str, output_path: str):
    output_with_format = os.path.join(output_path + '.json')
    # still need to filter just necessary information for json output
    # concatenated_data = {"interproscan-version": "6.0.0", 'results': data}
    with open(output_with_format, 'w') as json_file:
        json.dump(matches, json_file, indent=2)


def xml_output(matches: str, output_path: str):
    # A lot of changes! Need to be recreated!
    pass


def gff3_output(matches: str, output_path: str):
    pass


def write_results(matches_path: str, output_format: str, output_path: str):
    output_format = output_format.upper()
    with open(matches_path, 'r') as data:
        matches = json.load(data)
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
