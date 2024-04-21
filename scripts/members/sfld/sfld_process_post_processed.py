import argparse
import sys

from pathlib import Path


class SfldHit:
    def __init__(self):
        self.sequence = None  # query protein id
        self.sites = []   # store SiteMatch
        self.domains = []   # store domain InterPro model acc

    def add_feature(self, section: str, value: str):
        if section == "Domains:":
            # e.g.: SFLDF00315	0.000e+00	6.097e+02	0.300	1	443	609.600	1	447	1	448	0.000e+00	0.000e+00	0.990	0.300
            self.domains.append(value.split()[0])
        elif section == "Sites:":
            # e.g.: SFLDF00315 C129,C133,C136 Binds [4Fe-4S]-AdoMet cluster
            _site = SiteMatch()
            _site.query_ac = self.sequence
            _site.model_ac = value.split()[0]
            _site.site_residues = value.split()[1]
            _site.site_desc = " ".join(value.split()[2:])
            self.sites.append(_site)


class SiteMatch:
    def __init__(self):
        self.query_ac = None
        self.model_ac = None
        self.site_residues = None
        self.site_desc = None


def build_parser() -> argparse.ArgumentParser:
    """Build cmd-line argument parser"""
    parser = argparse.ArgumentParser(
        prog="process_post_processed_SFLD_hits",
        description="Filter results of HMMER search on SFLD HMMs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "sfld",
        type=Path,
        help="Path to post-processed SFLD results"
    )

    parser.add_argument(
        "-d", "--dom",
        type=Path,
        default=None,
        help="HMMER domtblout file",
    )

    parser.add_argument(
        "-O", "--hmmer_out",
        type=Path,
        default=None,
        help="HMMER output file",
    )

    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=sys.stdout,
        help="output file (otherwise STDOUT)",
    )

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    hits = parse_sfld(args.sfld)

    if args.dom:
        parse_dtbl(args.dom, hits)

    # if args.hmmer_out:
    #     print()


def parse_sfld(sfld: Path) -> list[SfldHit]:
    """Retrieve site and domain hits from post-processed SFLD output

    :param sfld: path to post-processed SFLD output

    :return: dict, keyed by query protein accession,
        values of SiteMatch instances and InterPro families
        (representing domain hits)
    """
    hits = {}  # {prot id : Hit}

    with open(sfld, "r") as fh:
        hit = SfldHit()
        section = None
        for line in fh:
            if line.strip() == "//":
                if hit.sequence:
                    try:
                        hits[hit.sequence]
                    except KeyError:
                        hits[hit.sequence] = hit
                hit = SfldHit()
            elif line.startswith("Sequence:"):
                prot_id = line.split()[-1].split("|")[-1]
                try:
                    hit = hits[prot_id]
                except KeyError:
                    hit.sequence = prot_id
            elif line.strip() in ("Domains:", "Sites:"):
                section = line.strip()
            else:
                hit.add_feature(section, line)

    try:
        hits[hit.sequence]
    except KeyError:
        hits[hit.sequence] = hit

    return hits


def parse_dtbl(dtbl: Path, hits: dict[str, SfldHit]):
    """Parse hmmer.dtbl output, removing hits not in 
    slfd-post-processed, and adding site data

    :param dtbl: Path to hmmer.dtbl file
    :param hits: dict of hits from sfld post-processed
    """
    processed_file = Path(dtbl.parent) / dtbl.name.replace(".dtbl", ".processed.dtbl")
    print(processed_file)
    with open(processed_file, "w") as out_fh:
        with open(dtbl, "r") as in_fh:
            for line in in_fh:
                if line.startswith("#"):
                    out_fh.write(line)
                else:
                    # check if domain hit in sfld processed
                    prot_id = line.split()[0]
                    model_ac = line.split()[4]

                    try:
                        if model_ac in hits[prot_id].domains:
                            out_fh.write(line)
                    except KeyError:
                        continue  # not in sfld-processed file

        # add site annotations
        out_fh.write("#\n# Sites\n# target name\taccession\tsite residue\tsite description\n")
        out_fh.write("#------------------- ---------- ------------ ----------------\n")
        for prot_id in hits:
            for _site in hits[prot_id].sites:
                out_fh.write(
                    f"{prot_id}\t{_site.model_ac}\t{_site.site_residues}\t{_site.site_desc}\n"
                )


if __name__ == "__main__":
    main()
