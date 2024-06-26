"""Parse output from pfsearchV3 wrapper run_prosite.py and
filter out blacklisted models"""
import sys
import json

from dataclasses import dataclass


@dataclass
class PrositeHit:
    profile: str = None
    profile_name: str = None
    motif_start: int = None
    motif_end: int = None
    sequence_id: str = None  # query prot seq id
    start: int = None  # start in prot seq
    end: int = None  # end in prot seq
    raw_score: float = None
    norm_score: float = None
    symbol: str = None  # always seems to be '+' - don't know why
    alignment: str = None

    @classmethod
    def from_list(cls, value: list[str]) -> 'PrositeHit':
        profile, profile_name = value[0].split("|")
        return cls(
            profile=profile,
            profile_name=profile_name,
            motif_start=int(value[1]),
            motif_end=int(value[2]),
            sequence_id=value[3],
            start=int(value[4]),
            end=int(value[5]),
            raw_score=float(value[6]),
            norm_score=float(value[7]),
            symbol=value[8],
            alignment=value[9]
        )


def _get_blacklist(blacklist_path: str) -> list[str]:
    """Retrieve the blacklisted profile IDs"""
    with open(blacklist_path, 'r') as fh:
        blacklist = fh.read().splitlines()
    return blacklist


def parse_and_filter(pfsearch_out: str, blacklist: list[str], release: str) -> dict:
    """Parse and filter the pfsearchV3 hits into the internal IPS6 JSON

    :param pfsearch_out: str repr of path to the output file from run_pfsearch.py
    :param blacklist: list of blacklisted profiles
    :param release: str of PROSITE Patterns release

    ips6_matches is keyed by the query protein sequence ID, and valued by
        a dict storing all hits for that seq. This dict is keyed by 
        signature accessions, and stores all matches against the signature/
        profile as "locations".
    """
    ips6_matches = {}
    with open(pfsearch_out, 'r') as fh:
        for line in fh:
            if line.strip():
                line_data = line.split()

                if line_data[0].split("|")[0] in blacklist:
                    continue

                hit = PrositeHit.from_list(line_data)

                if hit.sequence_id not in ips6_matches:
                    ips6_matches[hit.sequence_id] = {}

                if hit.profile not in ips6_matches[hit.sequence_id]:
                    ips6_matches[hit.sequence_id][hit.profile] = {
                        "accession": hit.profile,
                        "name": hit.profile_name,
                        "description": "",
                        "member_db": "PROSITE Profiles",
                        "version": release,
                        "model-ac": hit.profile,
                        "locations": []
                    }

                ips6_matches[hit.sequence_id][hit.profile]["locations"].append(
                    {
                        "start": hit.start,
                        "end": hit.end,
                        "representative": "",
                        "score": hit.norm_score,
                        "alignment": hit.alignment,
                        "location-fragments": [{
                            "start": hit.start,
                            "end": hit.end,
                            "dc-status": "CONTINUOUS"
                        }]
                    }
                )

    return ips6_matches


def main():
    pfsearch_out = sys.argv[1]  # output file from PFSEARCH_RUNNER module
    output_file = sys.argv[2]  # str rep of path for internal IPS6 JSON
    blacklist_path = sys.argv[3]  # str rep of path to file listing blacklisted models
    release = pfsearch_out.split("/")[-1].split("._.")[0].strip()

    blacklist = _get_blacklist(blacklist_path)

    ips6_matches = parse_and_filter(pfsearch_out, blacklist, release)

    with open(output_file, "w") as fh:
        json.dump(ips6_matches, fh, indent=2)


if __name__ == "__main__":
    main()
