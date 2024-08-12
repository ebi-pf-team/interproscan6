import json
import re
import sys

"""Matches up the corresponding IPS6 JSON file with the cath_resolve 
output file, and then runs filter_ips6_hits.py for each pair
of matches output files"""


CATH_PATTERN = re.compile(r"^.*\._\.(\d+\.\d+\.\d+\.\d+)\.cath\.resolved\.out$")
JSON_PATTERN = re.compile(r"^hmmer_parsed_.*\._\.(\d+\.\d+\.\d+\.\d+)\.json$")
RELEASE_PATTERN = re.compile(r"^\d+\.\d+\.\d+$")


def main():
    """
    Args include:
    All hmmer.out files from HMMER_PARSER
    All output files from cath-resolve hits
    Ending with the release number of FunFam from members.config
    """
    release = None
    files = {}  # keyed by cath superfamily
    for input_arg in sys.argv[1:]:
        _file = CATH_PATTERN.match(input_arg)
        if _file:
            cath_superfam = _file.group(1)
            if cath_superfam not in files:
                files[cath_superfam] = {}
            files[cath_superfam]["cath.resolve"] = input_arg
            continue

        _file = JSON_PATTERN.match(input_arg)
        if _file:
            cath_superfam = _file.group(1)
            if cath_superfam not in files:
                files[cath_superfam] = {}
            files[cath_superfam]["ips6.json"] = input_arg
            continue

        if RELEASE_PATTERN.match(input_arg):
            release = input_arg
            continue

        print(f"Did not recognise this input arg {input_arg}")
        sys.exit(1)

    for cath_superfam, file_info in files.items():
        cath_out = parse_cath(file_info["cath.resolve"])
        processed_ips6 = filter_matches(
            file_info["ips6.json"],
            cath_out,
            release
        )

        with open(f"{file_info['ips6.json']}.processed.json", "w") as fh:
            json.dump(processed_ips6, fh, indent=2)


class FunfamHit:
    def __init__(self):
        self.sequence_id = None
        self.domains = {}

    def add_domain(self, value: str):
        value = value.split()
        # carried over from i5:
        # treat each domain range in resolved hits as a separate domain
        for domain_range in value[4].split(","):
            match_id = value[1]
            domain = DomainHit()
            domain.signature_acc = match_id
            domain.score = value[2]
            domain.evalue = value[-1]
            domain.boundaries_start = domain_range.split("-")[0]
            domain.boundaries_end = domain_range.split("-")[-1]
            domain.resolved = value[4]

            if match_id not in self.domains:
                self.domains[match_id] = [domain]
            else:
                self.domains[match_id].append(domain)


class DomainHit:
    def __init__(self):
        self.signature_acc = None
        self.score = None  # bit score
        self.evalue = None  # indp-evalue
        self.boundaries_start = None  # resolved region
        self.boundaries_end = None  # resolved region
        self.resolved = None


def parse_cath(cath_out: str) -> dict[str, FunfamHit]:
    """Parse cath_out file into a dictionary

    :param cath_out: str repr of path to add_cath_superfamilies.py output file
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


def filter_matches(
    ips6: str,
    funfam_matches: dict[str, FunfamHit],
    release: str
) -> tuple[dict, set]:
    """Parse the IPS6 JSON file, filtering hits to only retains
    those that passed the Gene3D post-processing.

    :param ips6: str repr of path to internal IPS6 JSON file containing parsed hits from HMMER.out file
    :param funfam_matches: dict of FunfamHits, representing hits in the 
        add_cath_superfamilies.py output file
    :param release: FunFam release version

    Return processed IPS6 dict and a list of all cath superfamilies where hits were generated
    """
    processed_ips6 = {}
    with open(ips6, "r") as fh:
        ips6_data = json.load(fh)

    for protein_id in ips6_data:
        if protein_id not in funfam_matches:
            continue

        for signature_acc in ips6_data[protein_id]:  # e.g. 3.40.50.1170-FF-000001
            if signature_acc not in funfam_matches[protein_id].domains:
                continue

            # retrieve the relevant domain hit
            ips6_location = None  # from the IPS6 data
            funfam_domain = None  # from the cath-superfamilies output
            for location in ips6_data[protein_id][signature_acc]["locations"]:  # list of locations
                for domain in funfam_matches[protein_id].domains[signature_acc]:
                    # Iterate the list of DomainHit instances
                    if str(domain.evalue) == str(location["evalue"]) \
                        and float(domain.score) == float(location["score"]) \
                        and int(domain.boundaries_start) == int(location["envelopeStart"]) \
                        and int(domain.boundaries_end) == int(location["envelopeEnd"]):
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


if __name__ == "__main__":
    main()
