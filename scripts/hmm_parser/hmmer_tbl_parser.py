import argparse
import json

COMMENT_LINE = "#"


def parse(tbl_file, appl):
    tbl_info = []
    with open(tbl_file, "r") as f:
        for line in f.readlines():
            if line.startswith(COMMENT_LINE):
                pass
            else:
                tbl_info.append(line)


def main():
    parser = argparse.ArgumentParser(
        description="hmmer parser"
    )
    parser.add_argument(
        "-preproc", "--preproc_tbl", type=str, help="tbl file result of hmmer preproc")
    parser.add_argument(
        "-appl", "--application", type=str, help="name of member database")
    args = parser.parse_args()

    hmmer_parse_result = parse(args.preproc_tbl, args.application)
    print(json.dumps(hmmer_parse_result))


if __name__ == "__main__":
    main()
