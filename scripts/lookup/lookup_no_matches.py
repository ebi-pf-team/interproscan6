import argparse
import ast
import json


def main():
    parser = argparse.ArgumentParser(
        description="Parse to no match lookup"
    )
    parser.add_argument(
        "-checked", "--checked_lookup", type=str, help="dict with md5 lookup checked"
    )
    args = parser.parse_args()

    with open(args.checked_lookup, 'r') as md5_data:
        checked_data = md5_data.read()
    checked_info = ast.literal_eval(checked_data)
    no_matches = checked_info["no_matches"]

    if no_matches:
        print(no_matches)


if __name__ == "__main__":
    main()
