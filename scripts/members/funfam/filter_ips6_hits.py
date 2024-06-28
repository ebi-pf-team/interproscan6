import argparse
import json

from pathlib import Path


class FunfamHit:
    def __init__(self):
        self.sequence_id = None
        self.domains = {}

    def add_domain(self, value: str):
        match_id = value.split()[1]
        domain = DomainHit()
        domain.signature_acc = match_id
        domain.score = value.split()[2]
        domain.evalue = value.split()[-1]
        domain.boundaries_start = value.split()[3].split("-")[0]
        domain.boundaries_end = value.split()[3].split("-")[1]
        domain.resolved = value.split()[4]
        domain.aligned_regions = value.split()[-3]

        if match_id not in self.domains:
            self.domains[match_id] = [domain]
        else:
            self.domains[match_id].append(domain)


class DomainHit:
    def __init__(self):
        self.signature_acc = None
        self.score = None  # bit score
        self.evalue = None  # indp-evalue
        self.boundaries_start = None  # envelope boundaries
        self.boundaries_end = None
        self.resolved = None
        self.aligned_regions = None


def build_parser() -> argparse.ArgumentParser:
    """Build cmd-line argument parser"""
    parser = argparse.ArgumentParser(
        prog="fiter_ips6_hits",
        description="""
        Parse the output from the cath_resolve_hits,
        filtering and adding data to the hits in the interal IPS6 JSON structure
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "ips6",
        type=Path,
        help="Path to an internal IPS6 JSON file"
    )

    parser.add_argument(
        "cath_out",
        type=Path,
        help="Path to cath resolve output file"
    )

    parser.add_argument(
        "out_json",
        type=Path,
        help="Path to write the output JSON file"
    )

    parser.add_argument(
        "release",
        type=str,
        help="FunFam release version"
    )

    return parser


def parse_cath(cath_out: Path) -> dict[str, FunfamHit]:
    """Parse cath_out file into a dictionary

    :param cath_out: Path to add_cath_superfamilies.py output file
    """
    matches = {}  # prot seq id: domains

    with open(cath_out, "r") as fh:
        for line in fh.readlines():
            if line.startswith('#'):
                continue  # header row
            protein_id = line.split()[0]
            if protein_id not in matches:
                match = FunfamHit()
                match.sequence_id = protein_id
                matches[protein_id] = match
            matches[protein_id].add_domain(line)

    return matches


def filter_matches(ips6: Path, gene3d_matches: dict[str, FunfamHit], release: str) -> tuple[dict, set]:
    """Parse the IPS6 JSON file, filtering hits to only retains
    those that passed the Gene3D post-processing.

    :param ips6: path to internal IPS6 JSON file containing parsed hits from HMMER.out file
    :param gene3d_matches: dict of FunfamHits, representing hits in the 
        add_cath_superfamilies.py output file
    :param release: FunFam release version

    Return processed IPS6 dict and a list of all cath superfamilies where hits were generated
    """
    processed_ips6 = {}
    with open(ips6, "r") as fh:
        ips6_data = json.load(fh)

    for protein_id in ips6_data:
        if protein_id not in gene3d_matches:
            continue

        for signature_acc in ips6_data[protein_id]:  # e.g. 3.40.50.1170-FF-000001
            if signature_acc not in gene3d_matches[protein_id].domains:
                continue

            # retrieve the relevant domain hit
            ips6_location = None  # from the IPS6 data
            funfam_domain = None  # from the cath-superfamilies output
            for location in ips6_data[protein_id][signature_acc]["locations"]:  # list of locations
                for domain in gene3d_matches[protein_id].domains[signature_acc]:
                    # Iterate the list of DomainHit instances
                    if domain.evalue == location["evalue"] \
                        and domain.score == location["score"] \
                        and domain.boundaries_start == location["envelopeStart"] \
                        and domain.boundaries_end == location["envelopeEnd"]:
                        ips6_location = location  # dict of hmmer match data
                        funfam_domain = domain  # DomainHit instance

                if not ips6_location:
                    continue

                if protein_id not in processed_ips6:
                    processed_ips6[protein_id] = {}

                funfam_sig_acc = f"G3DSA:{signature_acc}".replace("-", ":")
                if funfam_sig_acc not in processed_ips6[protein_id]:
                    sig_info = ips6_data[protein_id][signature_acc]
                    sig_info["member_db"] = "funfam"
                    sig_info["version"] = release
                    sig_info["accession"] = funfam_sig_acc

                    # model ac is the domain id (minus the -... suffix)
                    sig_info["model-ac"] = funfam_sig_acc.replace(":", "-")

                    processed_ips6[protein_id][funfam_sig_acc] = sig_info
                    processed_ips6[protein_id][funfam_sig_acc]["locations"] = []
                    # start locations as empty as not all hits/locations in ips6 
                    # may have parsed the post-processing

                # add the location fragments (the 'aligned-regions') to the domain location data
                ips6_location["location-fragments"] = []
                if len(funfam_domain.resolved.split(",")) == 1:
                    ips6_location["location-fragments"].append({
                        "start": funfam_domain.resolved.split("-")[0],
                        "end": funfam_domain.resolved.split("-")[1],
                        "dc-status": "CONTINUOUS"
                    })
                else:
                    # first fragment has dc-status C_TERMINAL_DISC
                    # the last fragment, dc-status = N_TERMINAL_DISC
                    # all fragments in between = NC_TERMIANL_DISC
                    # Normally the fragments are listed in order (c-term to n-term), but best to check
                    all_frags = [(int(_.split("-")[0]), int(_.split("-")[1])) for _ in funfam_domain.resolved.split(",")]
                    all_frags = sorted(all_frags, key=lambda x: (x[0], x[1]))

                    ips6_location["location-fragments"].append({
                        "start": all_frags[0],
                        "end": all_frags[1],
                        "dc-status": "C_TERMINAL_DISC"
                    })

                    for fragment in all_frags[1:-1]:
                        ips6_location["location-fragments"].append({
                            "start": fragment[0],
                            "end": fragment[1],
                            "dc-status": "NC_TERMINAL_DISC"
                        })

                    ips6_location["location-fragments"].append({
                        "start": all_frags[0],
                        "end": all_frags[1],
                        "dc-status": "N_TERMINAL_DISC"
                    })

                processed_ips6[protein_id][funfam_sig_acc]["locations"].append(ips6_location)

    return processed_ips6


def main():
    parser = build_parser()
    args = parser.parse_args()

    matches = parse_cath(args.cath_out)
    processed_ips6 = filter_matches(args.ips6, matches, args.release)

    with open(args.out_json, "w") as fh:
        json.dump(processed_ips6, fh, indent=2)


if __name__ == "__main__":
    main()
