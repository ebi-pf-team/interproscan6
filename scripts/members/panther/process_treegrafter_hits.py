import argparse

from pathlib import Path


METADATA_LINES = (
    "# Program:",
    "# Version:",
    "# Pipeline mode:",
    "# Query file:",
    "# Target file:",
    "# Option settings:",
    "# Current dir:",
    "# Date:",
    "# [ok]"
)


def main():
    parser = build_parser()
    args = parser.parse_args()

    hits = parse_treegrafter(args)

    update_dtbl(args.hmmer, hits)


def build_parser() -> argparse.ArgumentParser:
    """Build cmd-line argument parser"""
    parser = argparse.ArgumentParser(
        prog="process_post_processed_panther_hits",
        description="Filter results of HMMER search on Pather HMMs by parsing the TreeGrafter output",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "treegrafter",
        type=Path,
        help="Path to TreeGrafter output file"
    )

    parser.add_argument(
        "hmmer",
        type=Path,
        help="Path to HMMER .dtbl output file"
    )

    return parser


def parse_treegrafter(treegrafter: Path) -> dict[str, list[str]]:
    """Parse TreeGrafter output into a dict.

    :param treegrafter: Path to treeGrafter output file

    Retrun dict, keyed by protein IDs and valued by list of Panther signature accessions
    """
    hits = {}
    with open(treegrafter, "r") as fh:
        for line in fh:
            if line.startswith("query_id"):
                # header row
                continue
            protein_id = line.split()[0]
            sig_acc = line.split()[1]
            try:
                hits[protein_id].append(sig_acc)
            except KeyError:
                hits[protein_id] = [sig_acc]

    return hits


def update_dtbl(dtbl: Path, hits: dict[str, list[str]]) -> None:
    """Parse hmmer.dtbl output, removing hits not in TreeGrafter output.

    :param dtbl: Path to hmmer.dtbl file
    :param hits: dict of hits from panther post-processed output
        (i.e. TreeGrafter output)
    """
    processed_file = Path(dtbl.parent) / dtbl.name.replace(".dtbl", ".dtbl.post_processed.dtbl")
    closing_lines = ["#\n"]
    with open(processed_file, "w") as out_fh:
        with open(dtbl, "r") as in_fh:
            for line in in_fh:
                if line.startswith("#"):
                    if line.startswith(METADATA_LINES):
                        closing_lines.append(line)
                    else:
                        out_fh.write(line)
                else:
                    # check if ther protein is in treegrafter output
                    prot_id = line.split()[0]
                    sig_acc = line.split()[3]
                    try:
                        sig_accessions = [_.split(":")[0] for _ in hits[prot_id]]
                        if sig_acc.split(".") in sig_accessions:
                            out_fh.write(line)
                    except KeyError:
                        # protein was not in the TreeGrafter output
                        continue

        closing_lines.append("#")
        for line in closing_lines:
            out_fh.write(line)


if __name__ == '__main__':
    main()