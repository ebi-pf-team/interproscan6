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
    subparsers = parser.add_subparsers()

    parser_pre = subparsers.add_parser(
        "prepare",
        help="convert the PAINT annotation file into per-family JSON files"
    )
    parser_pre.add_argument("annotation_file", type=Path,
                            help="PAINT annotation file "
                                 "(e.g. PAINT_Annotations_TOTAL.txt)")
    parser_pre.add_argument("output_dir", type=Path,
                            help="output directory for per-family JSON files")
    parser_pre.set_defaults(func=prepare)

    parser_run = subparsers.add_parser("run")
    parser_run.add_argument("-t", "--threads", type=int, default=1)
    parser_run.add_argument("jsonfile", type=Path)
    parser_run.add_argument("msfdir", type=Path)
    parser_run.set_defaults(func=run)

    args = parser.parse_args()
    try:
        func = args.func
    except AttributeError:
        parser.error("too few arguments")
    else:
        func(args)


def run(args):
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

        length = get_alignment_width(fasta_path)

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


def get_alignment_width(fasta_path: Path) -> int:
    width = 0
    in_first_sequence = False
    with fasta_path.open("rt") as fh:
        for line in map(str.rstrip, fh):
            if not line:
                continue
            if line.startswith(">"):
                if in_first_sequence:
                    break
                in_first_sequence = True
                continue
            if in_first_sequence:
                width += len(line)

    return width


def prepare(args):
    """Convert the PAINT annotation file into one JSON file per family."""
    args.output_dir.mkdir(parents=True, exist_ok=True)
    families = {}
    with args.annotation_file.open("rt") as fh:
        for i, line in enumerate(fh):
            fam_an_id, annotations, graft_point = line.rstrip().split("\t")
            fam_id, node_id = fam_an_id.split(":")
            fam = families.setdefault(fam_id, {})
            go_terms = []
            protein_class = subfam_id = None
            for annotation in re.split(r"\s+|;", annotations):
                if re.fullmatch(r"PTHR\d+:(SF\d+)", annotation):
                    subfam_id = annotation
                elif re.fullmatch(r"GO:\d{7}", annotation):
                    go_terms.append(annotation)
                elif re.fullmatch(r"PC\d{5}", annotation):
                    protein_class = annotation
            fam[node_id] = [
                subfam_id,
                ",".join(go_terms) if go_terms else None,
                protein_class,
                graft_point
            ]
    for fam_id, obj in families.items():
        with (args.output_dir / f"{fam_id}.json").open("wt") as fh:
            json.dump(obj, fh)


if __name__ == "__main__":
    main()