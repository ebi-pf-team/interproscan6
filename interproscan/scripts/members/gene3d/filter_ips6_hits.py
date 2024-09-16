import json
import sys

from pathlib import Path


class Gene3dHit:
    def __init__(self):
        self.sequence_id = None
        self.domains = {}

    def add_domain(self, value: str):
        """Where value is a result of line.split()"""
        domain_id = value[0]
        match_id = value[3].replace(f"_{domain_id}", "")
        domain = DomainHit()
        domain.signature_acc = domain_id
        domain.cath_superfamily = value[1]
        domain.match_id = match_id
        domain.score = value[4]
        domain.evalue = value[9]
        domain.boundaries_start = value[5].split("-")[0]
        domain.boundaries_end = value[5].split("-")[1]
        domain.resolved = value[6]
        domain.aligned_regions = value[7]

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
        self.start = None  # resolved.split("-")[0]
        self.end = None  # resolved.split("-")[1]
        self.boundaries_start = None  # envelope boundaries
        self.boundaries_end = None
        self.resolved = None
        self.aligned_regions = None


def parse_cath(cath_out: Path) -> dict[str, Gene3dHit]:
    """Parse cath_out file into a dictionary

    :param cath_out: Path to add_cath_superfamilies.py output file
    """
    matches = {}  # prot seq id: domains

    with open(cath_out, "r") as fh:
        for line in fh.readlines():
            if line.startswith('#') or not line.strip():
                continue  # header row or blank line
            protein_id = line.split()[2]
            if protein_id not in matches:
                match = Gene3dHit()
                match.sequence_id = protein_id
                matches[protein_id] = match
            matches[protein_id].add_domain(line.split())

    return matches


def filter_matches(
    ips6: Path,
    gene3d_matches: dict[str, Gene3dHit],
    funfam_dir: Path,
) -> tuple[dict, set]:
    """Parse the IPS6 JSON file, filtering hits to only retains
    those that passed the Gene3D post-processing.

    :param ips6: path to internal IPS6 JSON file containing parsed hits from HMMER.out file
    :param gene3d_matches: dict of Gene3dHits, representing hits in the
        add_cath_superfamilies.py output file
    :param funfam_dir: path to funfam data dir where all hmms are stored -
        so it can check if the hmm exists

    Return processed IPS6 dict and a list of all cath superfamilies where hits were generated
    """
    processed_ips6 = {}
    all_cath_superfamilies = set()

    with open(ips6, "r") as fh:
        ips6_data = json.load(fh)

    for protein_id in ips6_data:
        if protein_id not in gene3d_matches:
            continue

        for match_id in ips6_data[protein_id]:
            if match_id not in gene3d_matches[protein_id].domains:
                continue

            # retrieve the relevant domain hit
            ips6_location = None  # from the IPS6 data
            gene3d_domain = None  # from the cath-superfamilies output

            # locations is a list of dicts, each representing a location where the signature
            # matched a region of the query protein sequece
            for location in ips6_data[protein_id][match_id]["locations"]:
                for domain in gene3d_matches[protein_id].domains[match_id]:
                    if str(domain.evalue) == str(location["evalue"]) \
                        and float(domain.score) == float(location["score"]) \
                        and int(domain.boundaries_start) == int(location["start"]) \
                        and int(domain.boundaries_end) == int(location["end"]):
                        ips6_location = location
                        gene3d_domain = domain

                if not ips6_location:
                    continue

                if protein_id not in processed_ips6:
                    processed_ips6[protein_id] = {}

                cath_superfam = gene3d_domain.cath_superfamily

                # only take cath-superfam fowards for FunFam analysis if HMM exists
                hmm_path = funfam_dir / f"{cath_superfam.replace('.', '/')}.hmm"
                if hmm_path.is_file():
                    all_cath_superfamilies.add(cath_superfam)

                gene3d_sig_acc = f"G3DSA:{cath_superfam}"
                if gene3d_sig_acc not in processed_ips6[protein_id]:
                    # signature_acc is the domain id
                    # replace the domain id with the Cath superfamily
                    sig_info = ips6_data[protein_id][match_id]
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
    """CL input:
    0. Path to an internal IPS6 JSON file
    1. Path to cath_out file
    2. Path to write the output JSON file
    3. Path to write out a plain text file listing the CATH superfamilies
    4. Path to FunFam models dir, e.g. data/funfam/models
    """
    args = sys.argv[1:]
    ips6 = Path(args[0])
    cath_out = Path(args[1])
    out_json = Path(args[2])
    out_superfamilies = Path(args[3])
    funfam = Path(args[4])

    matches = parse_cath(cath_out)
    processed_ips6, superfamilies = filter_matches(ips6, matches, funfam)

    with open(out_json, "w") as fh:
        json.dump(processed_ips6, fh)
    with open(out_superfamilies, "w") as fh:
        for superfam in superfamilies:
            fh.write(f"{superfam}\n")


if __name__ == "__main__":
    main()
