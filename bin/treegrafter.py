#!/usr/bin/env python3


import argparse
import json
import re
import subprocess
from pathlib import Path

from Bio import Phylo
from Bio.Phylo import NewickIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("workdir", type=Path)
    parser.add_argument("msfdir", type=Path)
    args = parser.parse_args()

    for fasta in args.workdir.glob("*faa"):
        sequence_id, family_id = parse_header(fasta)

        jplace = run_epang(
            fasta,
            args.msfdir / f"{family_id}.AN.fasta",
            args.msfdir / f"{family_id}.bifurcate.newick",
            threads=args.threads
        )
        
        if jplace:
            tree = args.msfdir / f"{family_id}.newick"
            for query_id, node_id in parse_jplace(jplace, tree):
                print(sequence_id, family_id, node_id)


def run_epang(fastafile: Path, msafile: Path, treefile: Path, threads: int = 1) -> Path | None:
    proc = subprocess.run([
        "epa-ng",
        "-G", "0.05",
        "-m", "WAG",
        "-T", str(threads),
        "-t", str(treefile),
        "-s", str(msafile),
        "-q", str(fastafile),
        "--redo"
    ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    result_file = Path("epa_result.jplace")
    if proc.returncode == 0 and result_file.is_file():
        return result_file
    return None


def parse_jplace(jplacefile: Path, treefile: Path):
    with jplacefile.open("rt") as fh:
        results = json.load(fh)

    tree_string = results["tree"]
    matches = re.findall(r"AN(\d+):\d+\.\d+\{(\d+)\}", tree_string)

    an_label = {}
    for [an, r] in matches:
        an_label["AN" + an] = "R" + r
        an_label["R" + r] = "AN" + an

    newick_string = re.sub(r"(AN\d+)?\:\d+\.\d+{(\d+)}", r"R\g<2>",
                           tree_string)

    newick_string = re.sub(r"AN\d+", r"", newick_string)
    newick_string = re.sub(r"BI\d+", r"", newick_string)
    mytree = Phylo.read(NewickIO.StringIO(newick_string), "newick")
    
    for placement in results["placements"]:
        query_id = placement["n"][0]
        child_ids = []
        ter = []

        for maploc in placement["p"]:
            rloc = "R" + str(maploc[0])
            clade_obj = mytree.find_clades(rloc)

            node = next(clade_obj)
            ter.extend(node.get_terminals())
            comonancestor = mytree.common_ancestor(ter)

            for leaf in comonancestor.get_terminals():
                child_ids.append(an_label[leaf.name])

        newtree = Phylo.read(treefile, "newick")
        common_an = newtree.common_ancestor(child_ids)
        yield query_id, str(common_an) if common_an else "root"


def parse_header(fasta: Path) -> tuple[str, str]:
    with fasta.open("rt") as fh:
        header = next(fh).rstrip()[1:]
        seq_id, family_id = header.split("|")
        return seq_id, family_id


if __name__ == "__main__":
    main()
