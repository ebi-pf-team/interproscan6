import csv
import argparse

SERIAL_GROUP = "PROTEIN"


# Default for protein sequences are TSV, XML and GFF3, for nucleotide sequences GFF3 and XML.
def write_results(members_info: list[dict], output_format: list):
    if len(output_format) == 0:
        if SERIAL_GROUP == "PROTEIN":
            output_format = ["TSV", "XML", "GFF3"]
        else:
            output_format = ["XML", "GFF3"]
    else:
        output_format = list(map(lambda x: x.upper(), output_format))

    if 'TSV' in output_format:
        with open(f"lookup.tsv", 'w') as file:
            for entry in members_info:
                writer = csv.writer(file, delimiter='\t')
                writer.writerow(entry.values())
                writer.writerow("\n")
        # if 'XML' in output_format:
        #     root = ET.Element('results')
        #     for result in members_info:
        #         el = ET.SubElement(root, 'result')
        #         el.text = result
        #     tree = ET.ElementTree(root)
        #     tree.write(output_file)
        # if 'JSON' in output_format:
        #     json.dump(members_info, file)
        # if 'GFF3' in output_format:


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
