import argparse
from typing import Iterator, TextIO

from pyhmmer.easel import Alphabet, SequenceFile
from pyhmmer.hmmer import hmmsearch as pyhmmsearch
from pyhmmer.plan7 import HMMFile, TopHits


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="hmmsearch using pyhmmer")
    cutoff_group = parser.add_argument_group("CUTOFF (ONLY ONE ALLOWED)")
    cutoff_args = cutoff_group.add_mutually_exclusive_group()
    parser.add_argument(
        "-i", "--input", required=True, help="input fasta sequence file"
    )
    parser.add_argument(
        "-d", "--database", required=True, help="scan_sequences database file"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="output tab-delimited results file"
    )
    parser.add_argument(
        "-c",
        "--cpus",
        type=int,
        default=20,
        help="number of cpus to use in parallel (default: %(default)s)",
    )
    parser.add_argument(
        "--domain",
        action="store_true",
        help="use to switch to domain searching mode. This also switches the cutoff to the domain equivalent. (default: %(default)s)",
    )

    cutoff_args.add_argument(
        "--trusted",
        action="store_true",
        help="use trusted cutoffs (TC), equivalent to --cut_tc with hmmsearch",
    )
    cutoff_args.add_argument(
        "-T",
        "--bit-threshold",
        type=float,
        default=40,
        help="bitscore cutoff (default: %(default)s)",
    )

    return parser.parse_args()


def _write_full_length_results(
    fdst: TextIO,
    name: str,
    hmmname: str,
    evalue: float,
    bitscore: float,
    dom_evalue: float,
    dom_bitscore: float,
    dom_start: int,
    dom_end: int,
):
    msg = f"{name}\t{hmmname}\t{evalue:.2e}\t{bitscore:.2f}\n"
    fdst.write(msg)


def _write_domain_results(
    fdst: TextIO,
    name: str,
    hmmname: str,
    evalue: float,
    bitscore: float,
    dom_evalue: float,
    dom_bitscore: float,
    dom_start: int,
    dom_end: int,
):
    msg = f"{name}\t{hmmname}\t{dom_evalue:.2e}\t{dom_bitscore:.2f}\t{dom_start}\t{dom_end}\n"
    fdst.write(msg)


def write_output(allhits: Iterator[TopHits], fdst: TextIO, domain_search: bool):
    if domain_search:
        header = "query\tscan_sequences\tdom_e-value\tdom_bitscore\tstart\tend\n"
        writer = _write_domain_results
    else:
        header = "query\tscan_sequences\te-value\tbitscore\n"
        writer = _write_full_length_results

    fdst.write(header)
    for hits in allhits:
        for hit in hits.reported:
            name = hit.name.decode()
            hmmname = hits.query_name.decode()
            evalue = hit.evalue
            bitscore = hit.score
            dom_evalue = hit.best_domain.i_evalue
            dom_bitscore = hit.best_domain.score
            dom_start = hit.best_domain.env_from
            dom_end = hit.best_domain.env_to
            writer(
                fdst=fdst,
                name=name,
                hmmname=hmmname,
                evalue=evalue,
                bitscore=bitscore,
                dom_evalue=dom_evalue,
                dom_bitscore=dom_bitscore,
                dom_start=dom_start,
                dom_end=dom_end,
            )


def main():
    args = parse_args()
    seqfile = args.input
    hmmfile = args.database
    output = args.output
    cpus = args.cpus
    domain_search = args.domain
    if args.trusted:
        cutoff = {"bit_cutoffs": "trusted"}
    elif not domain_search:
        cutoff = {"T": args.bit_threshold}
    else:
        cutoff = {"domT": args.bit_threshold}

    alphabet = Alphabet.amino()
    with SequenceFile(seqfile, digital=True, alphabet=alphabet) as sfp, HMMFile(
        hmmfile
    ) as hfp:
        sequences = sfp.read_block()
        hmms = hfp.optimized_profiles() if hfp.is_pressed() else hfp

        allhits: Iterator[TopHits] = pyhmmsearch(hmms, sequences, cpus=cpus, **cutoff)

        with open(output, "w") as ofp:
            write_output(allhits, ofp, domain_search)


if __name__ == "__main__":
    main()
