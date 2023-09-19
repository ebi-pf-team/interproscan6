import argparse
import json


def build(parsed_sequence, parsed_tbl, parsed_out):
    agg_result = []
    with open(parsed_sequence, "r") as seq_file:
        with open(parsed_tbl, "r") as tbl_file:
            seq = json.load(seq_file)
            tbl = json.load(tbl_file)
            for match in tbl:
                match["sequence"] = seq[match["id"]]
                agg_result.append(match)
    return tbl

# def build(parsed_sequence, parsed_tbl, parsed_out):
#     agg_result = []
#
#     with open(parsed_sequence, "r") as seq_file:
#         seq = json.load(seq_file)
#     with open(parsed_tbl, "r") as tbl_file:
#         tbl = json.load(tbl_file)
#     with open(parsed_out, "r") as out_file:
#         out = json.load(out_file)
#
#     for match_tbl in tbl:
#         for match_out in out:
#             id = match_out["sequence_match"].keys()
#             if match_tbl["id"] == next(iter(id)):
#                 match_tbl["sequence"] = seq[match_tbl["id"]]
#                 match_tbl["out"] = match_out
#                 agg_result.append(match_tbl)
#     return tbl


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

    agg_result = build(args.parsed_sequence, args.parsed_tbl, args.parsed_out)
    print(json.dumps(agg_result))


if __name__ == "__main__":
    main()
