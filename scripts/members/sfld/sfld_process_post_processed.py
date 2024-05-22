import argparse
import json
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
        "ips6",
        type=Path,
        default=None,
        help="Path to internal IPS6 JSON structure (from HMMER_PARSER)",
    )

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    hits = parse_sfld(args.sfld)
    updated_ips6 = filter_matches_and_add_site(args.ips6, hits)
    print(json.dumps(updated_ips6, indent=2))


def parse_sfld(sfld: Path) -> dict[str, SfldHit]:
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


def filter_matches_and_add_site(ips6, hits):
    """Parse internal IPS6 JSON and remove matches not selected by SFLD post-processing
    and add in site annotations

    :param ips6: path to internal IPS6 JSON structure
    :param hits: dict of SfldHit instances.
        keyed by query protein accession, with
        values of SiteMatch instances and InterPro families
        (representing domain hits)
    """
    with open(ips6, "r") as fh:
        ips6_data = json.load(fh)

    for protein_id in ips6_data:
        if protein_id not in hits:
            # if there are only SFLD hits (i.e. no hits from other tools): delete the protein
            # else only delete the SFLD keys
            if all(signature_acc.startswith('sfld') for signature_acc in ips6_data[protein_id]):
                del ips6_data[protein_id]
                continue
            else:
                for signature_acc in ips6_data[protein_id]:
                    if signature_acc.startswith('sfld'):
                        del ips6_data[protein_id][signature_acc]

        else:
            for signature_acc in ips6_data[protein_id]:
                if signature_acc.upper().startswith('SFLD'):
                    if signature_acc not in hits[protein_id].domains:
                        # domain match did not parse the filtering of the post-processing
                        del ips6_data[protein_id][signature_acc]
                    else:
                        # add site data
                        # check each SfldSite instance as there can be multiple site hits
                        # for each signature accession in a protein sequence
                        for site in hits[protein_id].sites:
                            if site.model_ac == signature_acc:
                                if "sites" not in ips6_data[protein_id][signature_acc]["locations"]:
                                    
                                    site_positions = set()
                                    for position in site.site_residues.split(","):
                                        residues = position[1:].split("-")
                                        site_positions.update(map(int, residues))
                                    earliest_site, latest_site = int(min(site_positions)), int(max(site_positions))
                                   
                                    # find the relevant (domain) location
                                    for i, location in enumerate(ips6_data[protein_id][signature_acc]["locations"]):
                                        if int(location["start"]) <= earliest_site and int(location["end"]) >= latest_site:
                                            if "sites" not in ips6_data[protein_id][signature_acc]["locations"][i]:
                                                 ips6_data[protein_id][signature_acc]["locations"][i]["sites"] = []
                                            
                                            site_info = {
                                                "description": site.site_desc,
                                                "numLocations": len(site.site_residues.split(",")),
                                                "siteLocations": []
                                            }
                                            for site_location in site.site_residues.split(","):
                                                site_info['siteLocations'].append({
                                                    "start": site_location[0],
                                                    "end": site_location.split("-")[0][1:],
                                                    "residue": site_location[0],
                                                })
                                            ips6_data[protein_id][signature_acc]["locations"][i]["sites"].append(site_info)

                                            break

    return ips6_data


if __name__ == "__main__":
    main()
