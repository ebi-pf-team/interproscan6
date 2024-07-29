import json
import re
import sys


def matches2rules(matches_path: str, rules_hash: dict):
    with open(matches_path, "r") as matches:
        matches_info = json.load(matches)

    for sequence_id, domains in matches_info.items():
        for model_id, domain in domains.items():
            rule = rules_hash[model_id]

            domHits = []
            for location in domain["locations"]:
                sequence_id = sequence_id
                model_id = model_id
                hmm_from = location["hmmStart"]
                hmm_to = location["hmmEnd"]
                hmm_align = location["hmm_alignment"]
                seq_from = location["start"]
                seq_to = location["end"]
                seq_align = location["alignment"]

                map = map_hmm_to_seq(hmm_from, hmm_align, seq_align)

                rule_sites = []
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
                                    "end": pos['end'],
                                    "residue": pos['condition'],
                                    "start": pos['start']
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
                        "envelopeStart": 1  # always 1 in i5 but we have location["envelopeStart"]
                        "envelopeEnd": 2  # always 2 in i5 but we have location["envelopeEnd"]
                        'scope': rule['Scope'],
                    }
                    domHits.append(domHit)
            domain["locations"] = domHits

    return matches_info


def map_hmm_to_seq(hmm_pos, hmm, seq):
    """
    map base positions from alignment, from query HMM coords to (ungapped) target sequence coords
    arguments are hmm_from position, hmm_align_seq, query_align_seq
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
