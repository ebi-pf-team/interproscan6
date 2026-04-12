#!/usr/bin/env python3

import argparse
import json
import re
import subprocess
from pathlib import Path
from tempfile import mkstemp

from Bio import Phylo
from Bio.Phylo import NewickIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("jsonfile", type=Path)
    parser.add_argument("msfdir", type=Path)
    args = parser.parse_args()

    assert args.jsonfile.is_file()
    assert args.msfdir.is_dir()

    with args.jsonfile.open("rt") as fh:
        sequences = json.load(fh)

    for seq_id, matches in sequences.items():
        # Ensure we only have one family
        assert len(matches) == 1
        _, match = matches.popitem()

        # Ensure we only have one domain
        assert len(match["locations"]) == 1

        location = match["locations"][0]
        q_aln = location["queryAlignment"]
        t_aln = location["targetAlignment"]
        assert len(q_aln) == len(t_aln)

        # Get the expected length of the sequence
        family_id = match["modelAccession"]
        fasta_path = args.msfdir / f"{family_id}.AN.fasta"
        assert fasta_path.is_file()

        length = 0
        with fasta_path.open("rt") as fh:
            for line in map(str.rstrip, fh):
                if line[0] != ">":
                    length += len(line)

        # Init sequence, and pad N-terminal 
        sequence = "-" * (location["hmmStart"] - 1)

        # Build sequence
        t_aln = re.sub(r"[UO]", "X", t_aln, flags=re.I)
        for i, seq_char in enumerate(t_aln):
            hmm_char = q_aln[i]

            if hmm_char != ".":
                sequence += seq_char

        # Pad C-terminal
        assert len(sequence) <= length
        while len(sequence) < length:
            sequence += "-"

        fd, fasta_path = mkstemp(suffix=".faa")
        with open(fd, "wt") as fh:
            fh.write(f">{seq_id}\n")
            for i in range(0, len(sequence), 60):
                fh.write(f"{sequence[i:i+60]}\n")

        fasta_path = Path(fasta_path)
        jplace = run_epang(
            fasta_path,
            args.msfdir / f"{family_id}.AN.fasta",
            args.msfdir / f"{family_id}.bifurcate.newick",
            threads=args.threads
        )

        fasta_path.unlink()
        
        if jplace:
            tree = args.msfdir / f"{family_id}.newick"
            for query_id, node_id in parse_jplace(jplace, tree):
                print(seq_id, family_id, node_id)


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


if __name__ == "__main__":
    main()
