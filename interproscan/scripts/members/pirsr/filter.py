import json
import re
import sys


def matches2rules(matches_path: str, rules_hash: dict):
    with open(matches_path, "r") as matches:
        matches_info = json.load(matches)

    for sequence_id, domains in matches_info.items():
        for model_id, domain in domains.items():
            try:
                rule = rules_hash[model_id]
            except KeyError:
                rule = None

            domHits = []
            sorted_locations = sorted(domain["locations"], key=lambda x: (x["evalue"], -x["score"]))
            for location in sorted_locations:
                hmm_from = location["hmmStart"]
                hmm_to = location["hmmEnd"]
                hmm_align = location["hmm_alignment"]
                seq_from = location["start"]
                seq_to = location["end"]
                seq_align = location["alignment"]

                map = map_hmm_to_seq(hmm_from, hmm_align, seq_align)

                rule_sites = []
                if rule:
                    for grp, positions in rule['Groups'].items():
                        pass_count = 0
                        positions_parsed = []
                        for pos_num, pos in enumerate(positions):
                            condition = re.sub(r'[-()]', lambda x: {'-': '', '(': '{', ')': '}'}[x.group()],
                                               pos['condition'])
                            condition = condition.replace('x', '.')
                            query_seq = seq_align.replace('-', '')

                            if pos['hmmStart'] < len(map) and pos['hmmEnd'] < len(map):
                                target_seq = query_seq[map[pos['hmmStart']]:map[pos['hmmEnd']] + 1]
                            else:
                                target_seq = ''

                            if re.fullmatch(condition, target_seq):
                                pass_count += 1
                                if pos['start'] == 'Nter':
                                    pos['start'] = seq_from
                                if pos['end'] == 'Cter':
                                    pos['end'] = seq_to

                            positions_parsed.append({
                                'description': pos['desc'],
                                "group": int(pos['group']),
                                "hmmEnd": pos['hmmEnd'],
                                "hmmStart": pos['hmmStart'],
                                "label": pos['label'],
                                "numLocations": 1,  # always 1 on i5 (change to len(positions)?)
                                "siteLocations": [
                                    {
                                        "end": pos['end'],  # different value on i5
                                        "residue": pos['condition'],
                                        "start": pos['start']  # different value on i5
                                    }
                                ]
                            })

                        if pass_count == len(positions):
                            rule_sites.extend(positions_parsed)

                if rule_sites:
                    domHit = {
                        'score': location["score"],
                        'evalue': location["evalue"],
                        'hmmStart': hmm_from,
                        'hmmEnd': hmm_to,
                        'hmmAlign': hmm_align,
                        'start': seq_from,
                        'end': seq_to,
                        'alignment': seq_align,
                        'sites': rule_sites,
                        "representative": '',
                        "hmmLength": location["hmmLength"],
                        "envelopeStart": location["envelopeStart"],
                        "envelopeEnd": location["envelopeEnd"]
                    }
                    domHits.append(domHit)
            domain["locations"] = domHits
            domain["score"] = sorted_locations[0]["score"]
            domain["evalue"] = sorted_locations[0]["evalue"]

    return matches_info


def map_hmm_to_seq(hmm_pos, hmm, seq):
    """
    map base positions from alignment, from query HMM coords to (ungapped) target sequence coords
   :param hmm_pos: Is hmm_from from hmmer.out file
   :param hmm: Is hmm_align from hmmer.out file
   :param seq: Is seq alignment from hmmer.out file
    """
    seq_pos = 0
    map = [0]

    for i in range(0, hmm_pos):
        map[hmm_pos:] = [-1]

    for i in range(0, len(hmm)):
        map[hmm_pos:] = [seq_pos]

        if hmm[i:i + 1] != '.':
            hmm_pos += 1
        if seq[i:i + 1] != '-':
            seq_pos += 1

    return map


if __name__ == '__main__':
    args = sys.argv[1:]

    matches = args[0]
    rules_path = args[1]

    with open(rules_path) as rulesfile:
        rules_hash = json.load(rulesfile)

    result = matches2rules(matches, rules_hash)
    print(json.dumps(result, indent=4))
