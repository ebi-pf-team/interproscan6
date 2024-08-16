import json
import sys

from pathlib import Path


class HamapHit:
    def __init__(self):
        self.sequence_id = None
        self.domains = {}

    def add_domain(self, value: str):
        domain = DomainHit()
        sig_acc = value.split()[0].split("|")[0]
        domain.signature_acc = sig_acc
        domain.sig_name = value.split()[0].split("|")[1]
        domain.motif_start = value.split()[1]
        domain.motif_end = value.split()[2]
        domain.query_seq_id = value.split()[3]
        domain.start = value.split()[4]
        domain.end = value.split()[5]
        domain.raw_score = value.split()[6]
        domain.norm_score = value.split()[7]
        domain.symbol = value.split()[8]
        domain.seq = value.split()[9]

        if sig_acc not in self.domains:
            self.domains[sig_acc] = [domain]
        else:
            self.domains[sig_acc].append(domain)


class DomainHit:
    def __init__(self):
        self.signature_acc = None
        self.sig_name = None
        self.motif_start = None
        self.motif_end = None
        self.query_seq_id = None
        self.start = None
        self.end = None
        self.raw_score = None
        self.norm_score = None
        self.symbol = None  # do not know this represents, often '+'
        self.seq = None
        self.match_id = None


def parse_pf_out(pf_output: Path) -> dict[str, HamapHit]:
    """Parse pf_output file into a dictionary

    :param pf_output: Path to pfsearch_wrapper.py output file
    """
    matches = {}  # prot seq id: domains

    with open(pf_output, "r") as fh:
        for line in fh.readlines():
            if not line.strip():
                continue
            protein_id = line.split()[3]
            if protein_id not in matches:
                match = HamapHit()
                match.sequence_id = protein_id
                matches[protein_id] = match
            matches[protein_id].add_domain(line)

    return matches


def filter_matches(ips6: Path, hamap_matches: dict[str, HamapHit]) -> tuple[dict, set]:
    """Parse the IPS6 JSON file, filtering hits to only retains
    those that passed the Gene3D post-processing.

    :param ips6: path to internal IPS6 JSON file containing parsed hits from HMMER.out file
    :param hamap_matches: dict of HamapHits, representing hits in the
        pfsearch_wrapper.py output file

    Return processed IPS6 dict and a list of all cath superfamilies where hits were generated
    """
    processed_ips6 = {}

    with open(ips6, "r") as fh:
        ips6_data = json.load(fh)

    for protein_id in ips6_data:
        if protein_id not in hamap_matches:
            continue

        for signature_acc in ips6_data[protein_id]:
            if signature_acc not in hamap_matches[protein_id].domains:
                continue
            domain_matches = hamap_matches[protein_id].domains[signature_acc]

            hmmer_match = ips6_data[protein_id][signature_acc]

            if protein_id not in processed_ips6:
                processed_ips6[protein_id] = {}

            if signature_acc not in processed_ips6[protein_id]:
                processed_ips6[protein_id][signature_acc] = {}
                processed_ips6[protein_id][signature_acc]["accession"] = hmmer_match["accession"]
                processed_ips6[protein_id][signature_acc]["name"] = hmmer_match["name"]
                processed_ips6[protein_id][signature_acc]["description"] = hmmer_match["description"]
                processed_ips6[protein_id][signature_acc]["member_db"] = hmmer_match["member_db"]
                processed_ips6[protein_id][signature_acc]["model-ac"] = hmmer_match["model-ac"]
                processed_ips6[protein_id][signature_acc]["locations"] = []

            for domain in domain_matches:
                processed_ips6[protein_id][signature_acc]["name"] = domain.sig_name
                processed_ips6[protein_id][signature_acc]["locations"].append({
                    "start": domain.start,
                    "end": domain.end,
                    "representative": "false",
                    "score": domain.norm_score,
                    "alignment": domain.seq,
                    "location-fragments": [{
                        "start": domain.start,
                        "end": domain.end,
                        "dc-status": "CONTINUOUS",
                    }]
                })

    return processed_ips6


def main():
    """CL input:
    0. Internal IPS6 JSON
    1. pfsearch_wrapper.py output file
    2. Str repr of path for final output JSON
    """
    args = sys.argv[1:]
    ips6_json = args[0]
    pf_output = args[1]
    output = args[2]

    matches = parse_pf_out(pf_output)
    processed_ips6 = filter_matches(ips6_json, matches)

    with open(output, "w") as fh:
        json.dump(processed_ips6, fh)


if __name__ == "__main__":
    main()
