import argparse
import json

from pathlib import Path


class Gene3dHit:
    def __init__(self):
        self.sequence_id = None
        self.domains = {}

    def add_domain(self, value: str):
        match_id = value.split()[3]
        domain = DomainHit()
        domain.signature_acc = value.split()[0]
        domain.cath_superfamily = value.split()[1]
        domain.match_id = match_id
        domain.score = value.split()[4]
        domain.evalue = value.split()[-1]
        domain.boundaries_start = value.split()[-5].split("-")[0]
        domain.boundaries_end = value.split()[-5].split("-")[1]
        domain.resolved = value.split()[-4]
        domain.aligned_regions = value.split()[-3]

        if match_id not in self.domains:
            self.domains[match_id] = [domain]
        else:
            self.domains[match_id].append(domain)


class DomainHit:
    def __init__(self):
        self.signature_acc = None
        self.cath_superfamily = None
        self.match_id = None
        self.score = None  # bit score
        self.evalue = None  # indp-evalue
        self.boundaries_start = None  # envelope boundaries
        self.boundaries_end = None
        self.resolved = None
        self.aligned_regions = None


def build_parser() -> argparse.ArgumentParser:
    """Build cmd-line argument parser"""
    parser = argparse.ArgumentParser(
        prog="gene3d_and_funfam_match_parser",
        description="Parse the output from the add_cath_superfamilies.py into the interal IPS6 JSON structure",
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
        help="Path to cath_out file"
    )

    parser.add_argument(
        "out_json",
        type=Path,
        help="Path to write the output JSON file"
    )

    parser.add_argument(
        "out_superfamilies",
        type=Path,
        help="Path to write out a plain text file listing the CATH superfamilies"
    )

    return parser


def parse_cath(cath_out: Path) -> dict[str, Gene3dHit]:
    """Parse cath_out file into a dictionary

    :param cath_out: Path to add_cath_superfamilies.py output file
    """
    matches = {}  # prot seq id: domains

    with open(cath_out, "r") as fh:
        for line in fh.readlines():
            if line.startswith('#'):
                continue  # header row
            protein_id = line.split()[2]
            if protein_id not in matches:
                match = Gene3dHit()
                match.sequence_id = protein_id
                matches[protein_id] = match
            matches[protein_id].add_domain(line)

    return matches


def filter_matches(ips6: Path, gene3d_matches: dict[str, Gene3dHit]) -> tuple[dict, set]:
    """Parse the IPS6 JSON file, filtering hits to only retains
    those that passed the Gene3D post-processing.

    :param ips6: path to internal IPS6 JSON file containing parsed hits from HMMER.out file
    :param gene3d_matches: dict of Gene3dHits, representing hits in the 
        add_cath_superfamilies.py output file

    Return processed IPS6 dict and a list of all cath superfamilies where hits were generated
    """
    processed_ips6 = {}
    all_cath_superfamilies = set()
    with open(ips6, "r") as fh:
        ips6_data = json.load(fh)
    for protein_id in ips6_data:
        if protein_id not in gene3d_matches:
            continue

        for signature_acc in ips6_data[protein_id]:
            if signature_acc not in gene3d_matches[protein_id].domains:
                continue

            # retrieve the relevant domain hit
            ips6_location = None  # from the IPS6 data
            gene3d_domain = None  # from the cath-superfamilies output
            for location in ips6_data[protein_id][signature_acc]["locations"]:
                for domain in gene3d_matches[protein_id].domains[signature_acc]:
                    if domain.evalue == location["evalue"] \
                        and domain.score == location["score"] \
                        and domain.boundaries_start == location["envelopeStart"] \
                        and domain.boundaries_end == location["envelopeEnd"]:
                        ips6_location = location
                        gene3d_domain = domain

                if not ips6_location:
                    continue

                if protein_id not in processed_ips6:
                    processed_ips6[protein_id] = {}

                cath_superfam = gene3d_domain.cath_superfamily
                all_cath_superfamilies.add(cath_superfam)
                gene3d_sig_acc = f"G3DSA:{cath_superfam}"
                if gene3d_sig_acc not in processed_ips6[protein_id]:
                    # signature_acc is the domain id
                    # replace the domain id with the Cath superfamily
                    sig_info = ips6_data[protein_id][signature_acc]
                    sig_info["accession"] = gene3d_sig_acc

                    # model ac is the domain id (minus the -... suffix)
                    sig_info["model-ac"] = gene3d_domain.signature_acc.split("-")[0]

                    processed_ips6[protein_id][gene3d_sig_acc] = sig_info
                    processed_ips6[protein_id][gene3d_sig_acc]["locations"] = []
                    # start locations as empty as not all hits/locations in ips6 
                    # may have parsed the post-processing

                # add the location fragments (the 'aligned-regions') to the domain location data
                ips6_location["location-fragments"] = []
                if len(gene3d_domain.resolved.split(",")) == 1:
                    ips6_location["location-fragments"].append({
                        "start": gene3d_domain.resolved.split("-")[0],
                        "end": gene3d_domain.resolved.split("-")[1],
                        "dc-status": "CONTINUOUS"
                    })
                else:
                    # first fragment has dc-status C_TERMINAL_DISC
                    # the last fragment, dc-status = N_TERMINAL_DISC
                    # all fragments in between = NC_TERMIANL_DISC
                    # Normally the fragments are listed in order (c-term to n-term), but best to check
                    all_frags = [(int(_.split("-")[0]), int(_.split("-")[1])) for _ in gene3d_domain.resolved.split(",")]
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

                processed_ips6[protein_id][gene3d_sig_acc]["locations"].append(ips6_location)

    return processed_ips6, all_cath_superfamilies


def main():
    parser = build_parser()
    args = parser.parse_args()

    matches = parse_cath(args.cath_out)
    processed_ips6, superfamilies = filter_matches(args.ips6, matches)

    with open(args.out_json, "w") as fh:
        json.dump(processed_ips6, fh, indent=2)
    with open(args.out_superfamilies, "w") as fh:
        for superfam in superfamilies:
            fh.write(f"{superfam}\n")


if __name__ == "__main__":
    main()
