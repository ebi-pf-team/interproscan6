import argparse
import json


def build(parsed_sequence, parsed_tbl):
    agg_result = []
    with open(parsed_sequence, "r") as seq_file:
        with open(parsed_tbl, "r") as tbl_file:
            seq = json.load(seq_file)
            tbl = json.load(tbl_file)
            for match in tbl:
                for seq_info in seq:
                    if seq_info["id"] == match["id"]:
                        match["sequence"] = seq_info
                        agg_result.append(match)
    return tbl


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

    agg_result = build(args.parsed_sequence, args.parsed_tbl)
    print(json.dumps(agg_result))


if __name__ == "__main__":
    main()
