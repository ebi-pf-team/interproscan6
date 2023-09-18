import argparse
import json


def build(parsed_sequence, parsed_tbl, parsed_out):
    parsed_info = []
    with open(parsed_sequence, "r") as seq_file:
        with open(parsed_tbl, "r") as tbl_file:
            for seq in seq_file.readlines():
                parsed_info.append(seq)
                for tbl in tbl_file.readlines():
                    if tbl == "null":
                        break
                    else:
                        parsed_info.append(tbl)
    return parsed_info


def main():
    parser = argparse.ArgumentParser(
        description="hmmer parser"
    )
    parser.add_argument(
        "-seq", "--parsed_sequence", type=str, help="sequences parsed")
    parser.add_argument(
        "-tbl", "--parsed_tbl", type=str, help="tbl file result parsed")
    parser.add_argument(
        "-out", "--parsed_out", type=str, help="out file result parsed")
    args = parser.parse_args()

    parsed_result = build(args.parsed_sequence, args.parsed_tbl, args.parsed_out)
    print(parsed_result)


if __name__ == "__main__":
    main()
