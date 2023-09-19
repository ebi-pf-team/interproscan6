import argparse
import json

COMMENT_LINE = "#"


def parse(tbl_file, appl):
    tbl_output = []
    with open(tbl_file, "r") as f:
        for line in f.readlines():
            tbl_info = {}
            if line.startswith(COMMENT_LINE):
                pass
            else:
                info = line.split()
                if info[0] != "None":
                    tbl_info["id"] = info[0]
                    tbl_info["tbl"] = info
                    tbl_info["appl"] = appl
            if tbl_info:
                tbl_output.append(tbl_info)
    return tbl_output


def main():
    parser = argparse.ArgumentParser(
        description="hmmer parser"
    )
    parser.add_argument(
        "-preproc", "--preproc_tbl", type=str, help="tbl file result of hmmer preproc")
    parser.add_argument(
        "-appl", "--application", type=str, help="name of member database")
    args = parser.parse_args()

    tbl_parsed = parse(args.preproc_tbl, args.application)
    print(json.dumps(tbl_parsed))


if __name__ == "__main__":
    main()
