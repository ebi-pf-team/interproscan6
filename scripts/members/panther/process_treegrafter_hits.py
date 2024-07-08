import argparse
import json

from copy import copy

from pathlib import Path


class PantherHit:
    def __init__(self, sequence=None):
        self.sequence = sequence  # query protein id
        self.domains = []
        self.signatures = set()

    def add_feature(self, line: str):
        line_parts = line.split()
        self.signatures.add(line_parts[1].split(":")[0])
        self.domains.append(
            {
                "signature-acc:superfamily": line_parts[1],
                "score": line_parts[2],
                "evalue": line_parts[3],
                "dom_score": line_parts[4],
                "dom_evalue": line_parts[5],
                "hmm_start": line_parts[6],
                "hmm_end": line_parts[7],
                "ali_start": line_parts[8],
                "ali_end": line_parts[9],
                "env_start": line_parts[10],
                "env_end": line_parts[11],
                "node_id": line_parts[-1]
            }
        )


def main():
    parser = build_parser()
    args = parser.parse_args()

    hits = parse_treegrafter(args.treegrafter)
    updated_ips6 = update_ips6(args.ips6, hits, args.paint_annnotations)
    print(json.dumps(updated_ips6, indent=2))


def build_parser() -> argparse.ArgumentParser:
    """Build cmd-line argument parser"""
    parser = argparse.ArgumentParser(
        prog="process_post_processed_panther_hits",
        description="Filter results of HMMER search on Pather HMMs by parsing the TreeGrafter output",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "treegrafter",
        type=Path,
        help="Path to TreeGrafter output file"
    )

    parser.add_argument(
        "ips6",
        type=Path,
        help="Path to IPS6 JSON file containing parsed content from HMMER.out file"
    )

    parser.add_argument(
        "paint_annnotations",
        type=Path,
        help="Path to PAINT annotations file from the Panther release"
    )

    return parser


def parse_treegrafter(treegrafter: Path) -> dict[str, list[str]]:
    """Parse TreeGrafter output into a dict.

    Example:
    {'OUTROSEMLOOKUP': [('PTHR43780:SF2', 'AN233')]}

    :param treegrafter: Path to treeGrafter output file

    Retrun dict, keyed by protein IDs and valued by list of Panther signature accessions
    """
    hits = {}
    with open(treegrafter, "r") as fh:
        for line in fh:
            if line.startswith("query_id"):  # header row
                continue
            protein_id = line.split()[0]
            hits.setdefault(protein_id, PantherHit(protein_id)).add_feature(line)
    return hits


def update_ips6(ips6: Path, hits: dict[str, list[str]], paint_anno_path: Path) -> None:
    """Parse ips6 json containing hits from the hmmer.out file,
    removing hits not in TreeGrafter output.

    :param ips6: Path to IPS6 JSON file
    :param hits: dict of hits from panther post-processed output
        (i.e. TreeGrafter output)
        {protein_id: [(sig_acc, node_id)]}
    :param paint_anno_path: path to Panther release Paint Annotation file
    """
    with open(ips6, "r") as fh:
        ips6_data = json.load(fh)

    processed_ips6_data = {}

    for _protein_id in ips6_data:
        protein_id = _protein_id.replace(".", "_")
        if protein_id not in hits:
            # IPS6 will only contain Panther hits at this stage
            # so don't need to check if need to retain hits from other tools
            continue

        else:
            for signature_acc in ips6_data[_protein_id]:

                if signature_acc not in hits[protein_id].signatures:
                    # signature hit did not parse TreeGrafter post-processing
                    continue

                for panther_hit in hits[protein_id].domains:
                    if panther_hit["signature-acc:superfamily"].split(":")[0] == signature_acc:

                        if protein_id not in processed_ips6_data:
                            processed_ips6_data[_protein_id] = {}
                        if signature_acc not in processed_ips6_data[_protein_id]:
                            processed_ips6_data[_protein_id][signature_acc] = copy(ips6_data[_protein_id][signature_acc])
                            processed_ips6_data[_protein_id][signature_acc]["locations"] = []
                        
                        # Panther/Treegrafter only has one (the best) match per protein
                        for location in ips6_data[_protein_id][signature_acc]["locations"]:
                            if panther_hit["ali_start"] == location["start"] and \
                                panther_hit["ali_end"] == location["end"] and \
                                panther_hit["hmm_start"] == location["hmmStart"] and \
                                panther_hit["hmm_end"] == location["hmmEnd"] and \
                                panther_hit["env_start"] == location["envelopeStart"] and \
                                panther_hit["env_end"] == location["envelopeEnd"]:
                                
                                processed_ips6_data[_protein_id][signature_acc]["locations"] = [{
                                    "start": panther_hit["ali_start"],
                                    "end": panther_hit["ali_end"],
                                    "representative": location["representative"],
                                    "hmmStart": panther_hit["hmm_start"],
                                    "hmmEnd": panther_hit["hmm_end"],
                                    "hmmLength": location["hmmLength"],
                                    "rawHmmBounds": location["rawHmmBounds"],
                                    "hmmBounds": location["hmmBounds"],
                                    "evalue": panther_hit["evalue"],
                                    "score": panther_hit["score"],
                                    "envelopeStart": panther_hit["env_start"],
                                    "envelopeEnd": panther_hit["env_end"],
                                    "postProcessed": location["postProcessed"],
                                    "alignment": location["alignment"],
                                    "cigar_alignment": location["cigar_alignment"]
                                }]

                        processed_ips6_data[_protein_id][signature_acc]["node_id"] = panther_hit["node_id"]
                        processed_ips6_data[_protein_id][signature_acc]["accession"] = panther_hit["signature-acc:superfamily"]
                        processed_ips6_data[_protein_id][signature_acc]["model-ac"] = panther_hit["signature-acc:superfamily"]

    return processed_ips6_data


if __name__ == '__main__':
    main()
