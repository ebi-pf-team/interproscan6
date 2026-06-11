#!/usr/bin/env python3

import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("fastafile", type=Path)
    parser.add_argument("jsonfile", type=Path)
    parser.add_argument("funfamroot", type=Path)
    args = parser.parse_args()
    assert args.fastafile.is_file()
    assert args.jsonfile.is_file()
    assert args.funfamroot.is_dir()

    with args.jsonfile.open("rt") as fh:
        sequences = json.load(fh)

    seen = set()
    for matches in sequences.values():
        for match in matches.values():
            accession = match["signature"]["accession"]
            assert accession is not None
            assert accession.startswith("G3DSA:")
            cath_id = accession[6:]
            if cath_id not in seen:
                seen.add(cath_id)
                c, a, t, h = cath_id.split(".")
                hmmfile = (args.funfamroot / c / a / t / h).with_suffix(".hmm")
                if hmmfile.is_file():
                    print(
                        "hmmsearch",
                        "-Z",
                        "65245",
                        "--cut_tc",
                        "--cpu",
                        str(args.threads),
                        str(hmmfile),
                        str(args.fastafile),
                        sep=" "
                    )


if __name__ == "__main__":
    main()
