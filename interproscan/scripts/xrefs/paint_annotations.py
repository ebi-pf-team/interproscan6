import json
import sys
from pathlib import Path


def add_paint_annotations(matches_path: str, paint_anno_dir: Path) -> dict:
    """Retrieve PAINT Annotations for Panther hits.

    Retrieve for pre-calc and calculated matches because PAINT
    annotations are not retrieved from the Match Lookup Service.

    :param matches_path: str repr of path to the internal IPS6 json
    :param paint_anno_dir: Path repr of path to the PAINT annotation dir
    """
    with open(matches_path, "r") as fh:
        matches = json.load(fh)

    for seq_id, match_info in matches.items():
        for sig_acc, data in match_info.items():
            if data["member_db"].upper() == "PANTHER":
                anno_path = paint_anno_dir / f"{sig_acc}.json"
                with open(anno_path, 'r') as fh:
                    paint_annotations = json.load(fh)
                    node_data = paint_annotations[data["node_id"]]

                matches[seq_id][sig_acc]["proteinClass"] = node_data[2]
                matches[seq_id][sig_acc]["graftPoint"] = node_data[3]

    return matches


def main():
    """CL input:
    0. Str repr of the path to the internal IPS6 JSON file
    1. Str repr of the paint annotation dir
    2. Str repr of the path for the output file"""
    args = sys.argv[1:]

    matches = args[0]
    paint_anno_dir = Path(args[1])

    parsed_matches = add_paint_annotations(matches, paint_anno_dir)

    with open(args[2], "w") as fh:
        json.dump(parsed_matches, fh)


if __name__ == "__main__":
    main()
