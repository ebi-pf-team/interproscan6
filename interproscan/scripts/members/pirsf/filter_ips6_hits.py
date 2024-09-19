import json
import re
import sys

from pathlib import Path

"""The PIRFT.pl script from PIRSF was replaced with an in-house pirsf.pl 
script in 2019 in i5 when a bug fix was required. That perl script was 
replaced with this Python script in IPS6."""


CHILD_DAT_LINE = re.compile(r">PIRSF\d+\schild:\s(.*)$")


class DatEntry:
    def __init__(self):
        self.model_id = None
        self.name = None
        self.mean_l = None
        self.std_l = None
        self.min_s = None
        self.mean_s = None
        self.std_s = None
        self.blast = None
        self.children = []

    def add_model_id(self, value: str):
        self.model_id = value.split()[0].strip(">")
        line_match = CHILD_DAT_LINE.match(value)
        if line_match:
            self.children = line_match.group(1).split()

    def add_name(self, value: str):
        self.name = value

    def add_data_values(self, value: str):
        self.mean_l, self.std_l, self.min_s, self.mean_s, self.std_s = [float(_) for _ in value.split()]

    def add_blast(self, value: str):
        self.blast = 0 if value.split()[1] == "NO" else 1


class pirsfHit:
    def __init__(self):
        self.model_id = None
        self.score = None
        self.data = None
        self.children = []

    def add_model_data(self, model_id: str, match_data: dict):
        self.model_id = model_id
        self.score = match_data["score"]
        self.data = match_data

    def add_child(self, child_id: str):
        self.children.append(child_id)


def get_signature_lengths(dtbl_path: str) -> dict[str, int]:
    """Retrieve the tlen (signature length) from the dtbl output file from HMMscan
    because PIRSF uses the tlen as the HmmLength"""
    sig_lengths = {}  # signature acc: tlen
    with open(dtbl_path, "r") as reader:
        for line in reader:
            if line.startswith('#'):
                continue
            line_segments = line.split()
            sig_lengths[line_segments[1]] = int(line_segments[2])
    return sig_lengths


def load_dat(dat_path: str) -> tuple[dict, dict]:
    """Load data from pirsf.dat file

    :param dat_path: str repr of path to pirsf.dat file

    Return dict with all data entries
    and dict with {subfamily: parent family}
    """
    stage = ''
    entry = None
    dat_entries = {}  # model id: datEntry
    children = {}  # subfamily: parent family
    with open(dat_path, "r") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if entry:
                    dat_entries[entry.model_id] = entry
                    for child in entry.children:
                        children[child] = entry.model_id
                entry = DatEntry()
                entry.add_model_id(line)
                stage = 'GET_NAME'
            elif stage == 'GET_NAME':
                entry.add_name(line)
                stage = 'GET_VALUES'
            elif stage == 'GET_VALUES':
                entry.add_data_values(line)
                stage = 'GET_BLAST'
            elif stage == 'GET_BLAST':
                entry.add_blast(line)
                stage = 'GET_MODEL_ID'

    if entry:
        dat_entries[entry.model_id] = entry
        for child in entry.children:
            children[child] = entry.model_id

    return dat_entries, children


def get_location_data(match: dict) -> tuple[int, int, int, int, float]:
    """A PIRSF model often match multiple locations in a given query sequence.
    Combine those locations into one, taking the largest range possible across
    the match locations, and sum the scores of all the matches together.

    These values are used for deciding whether to keep or reject the
    signature-proteinSeq match.
    """
    seq_start, seq_end, hmm_start, hmm_end, score = [0]*5
    for location in match["locations"]:
        # Yes, the post-processing perl script for PIRSF uses the envelope start/end as
        # the seq start and end
        seq_start = int(location["start"])   # ali from
        seq_end = int(location["end"])   # ali to
        hmm_start = location["hmmStart"]   # hmm from
        hmm_end = location["hmmEnd"]   # hmm to
        score += float(location["score"])

        # update to match georgetown 2017 script
        if int(location["start"]) < seq_start and location["hmmStart"] < hmm_start:
            seq_start = location["start"]
            hmm_start = location["hmmStart"]
        if int(location["end"]) < seq_end and location["hmmEnd"] < hmm_end:
            seq_end = location["end"]
            hmm_end = location["hmmEnd"]

    return seq_start, seq_end, hmm_start, hmm_end, score


def get_best_match(filtered_models: dict[str, pirsfHit], sig_lengths: dict[str, int]) -> dict:
    """Out of the filtered models, select the best match and it's associated subfamily,
    and only retain this model/family.
    """
    processed_match = {}

    # sort the matches by their score, higest to lowest
    filtered_models = dict(
        sorted(filtered_models.items(), key=lambda item: item[1].score, reverse=True)
    )
    for pirsf_acc in filtered_models:
        # we iterate here to get to the first acc that is not a subfamily
        if pirsf_acc.startswith("PIRSF5"):  # ignore if a sub-family
            continue

        processed_match[pirsf_acc] = filtered_models[pirsf_acc].data
        # there is only one location as the PIRSF post-processing combines all 
        # locations into a single location
        # and PIRFT usese the tlen value from the dtbl file as the hmmLength
        processed_match[pirsf_acc]["locations"][0]["hmmLength"] = sig_lengths[pirsf_acc]

        # if there is a child subfamily that passes the earlier filtering
        # also keep the best child subfamily
        child_models = {}
        for pirsf_subfam in filtered_models[pirsf_acc].children:
            child_models[pirsf_subfam] = filtered_models[pirsf_subfam]
        # sort children by score
        child_models = dict(
            sorted(child_models.items(), key=lambda item: item[1].score, reverse=True)
        )
        for pirsf_subfam in child_models:
            processed_match[pirsf_subfam] = filtered_models[pirsf_subfam].data
            processed_match[pirsf_subfam]["locations"][0]["hmmLength"] = sig_lengths[pirsf_subfam]
            break  # we only want the best

        break  # we only want the best match!

    return processed_match


def filter_matches(ips6_json: Path, sig_lengths: dict, pirsf_dat: dict, children: dict) -> dict:
    """
    For each protein in matches, 
        Filter the models to retain only those that pass the threshold
        Then select the best of these filtered models to be kept

    :param ips6_json: path to internal IPS6 JSON file with hmmer hits
    :param sig_lengths: tlen values from the hmmscan dtbl file, signature length
    :param pirsf_dat: dict containing data from the pirsf.dat file
    :param children: dict containing {subfam: parent fam}
    """
    with open(ips6_json, "r") as fh:
        matches = json.load(fh)
    processed_matches = {}

    for protein_id in matches:
        filtered_models = {}  # models that passes the filtering criteria

        # filter the models
        for model_id in matches[protein_id]:
            seq_len = matches[protein_id][model_id]["qlen"]
            seq_start, seq_end, hmm_start, hmm_end, score = get_location_data(
                matches[protein_id][model_id]
            )

            # calculate ratios and deviations
            ovl = (seq_end - seq_start + 1) / seq_len
            ratio = (hmm_end - hmm_start + 1) / (seq_end - seq_start + 1)
            ld = abs(seq_len - pirsf_dat[model_id].mean_l)

            # check if fam (has children) or subfamily (does not have children)
            if model_id in children:
                # If a sub-family, process slightly differently. Only consider the score.
                if ratio > 0.67 and score >= pirsf_dat[model_id].min_s:

                    parent_model_id = children[model_id]  # parent = parent model id

                    # if we have a subfamily match we should consider the parent to also be a match
                    # but check if there is a match for the parent
                    # if the there is no match for the parent, drop the match for the subfam
                    try:
                        pirf_hit = pirsfHit()
                        pirf_hit.add_model_data(parent_model_id, matches[protein_id][parent_model_id])
                        pirf_hit.add_child(model_id)
                        filtered_models[parent_model_id] = pirf_hit

                    except KeyError:  # no match for parent so continue
                        continue

                    pirf_hit = pirsfHit()
                    pirf_hit.add_model_data(model_id, matches[protein_id][model_id])
                    filtered_models[model_id] = pirf_hit

            elif ratio > 0.67 and \
                ovl >= 0.8 and \
                    (score >= pirsf_dat[model_id].min_s) and \
                    (ld < 3.5 * pirsf_dat[model_id].std_l or ld < 50):#
                # everything passes the threshold of length, score and standard deviations of length
                pirf_hit = pirsfHit()
                pirf_hit.add_model_data(model_id, matches[protein_id][model_id])
                filtered_models[model_id] = pirf_hit

        # keep only the best family and associated subfamily for this protein seq
        if filtered_models:
            processed_matches[protein_id] = get_best_match(filtered_models, sig_lengths)

    return processed_matches


def main():
    """CL arguments:
    :args[0]: str repr of path to the internal IPS6 JSON container HMMER matches
    :args[1]: str repr of path to the HmmScan dtbl file
    :args[2]: str repr of path to the PIRSF.dat file
    :args[3]: str repr of the output path
    """
    args = sys.argv[1:]
    ips6 = args[0]
    dtbl = args[1]
    dat = args[2]
    out_json = args[3]

    sig_lenths = get_signature_lengths(dtbl)
    pirsf_dat, children = load_dat(dat)
    # children {subfam: parent}
    print(len(pirsf_dat), len(children))
    with open("children", "w") as fh:
        json.dump(children, fh)
    with open("pirsf_dat", "w") as fh:
        fh.write(str(pirsf_dat))
    print(pirsf_dat['PIRSF001220'])
    procossed_results = filter_matches(ips6, sig_lenths, pirsf_dat, children)

    with open(out_json, "w") as fh:
        json.dump(procossed_results, fh, indent=2)


if __name__ == "__main__":
    main()
