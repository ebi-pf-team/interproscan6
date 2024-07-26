#!/usr/bin/env python

import json
import re
import sys


def matches2rules(matches, rule, member_db, version):
    with open(matches, "r") as matches:
        matches_info = json.load(matches)

    for protein_id, domains in matches_info.items():
        sequence_id = protein_id
        for model_id, domain in domains.items():
            rule = rules_hash[model_id]
            for location in domain["locations"]:
                hmm_from = location["hmmStart"]
                hmm_to = location["hmmEnd"]
                hmm_align = location["hmmAlign"]
                seq_from = location["start"]
                seq_to = location["end"]
                hmm_align = location["alignment"]

                map = map_hmm_to_seq(hmm_from, hmm_align, seq_align)

                rule_sites = []

                for grp in rule['Groups'].keys():

                    # if model_id in matches[sequence_id]:
                    #     next

                    pass_count = 0

                    pos_num = -1
                    for pos in rule['Groups'][grp]:
                        pos_num += 1
                        condition = pos['condition']

                        condition = re.sub('-', '', condition)
                        condition = re.sub('\\(', '{', condition)
                        condition = re.sub('\\)', '}', condition)
                        condition = re.sub('x', '.', condition)

                        # query_seq = re.sub('-', '', domain[]["alignment"])

                        if pos['hmmStart'] < len(map) and pos['hmmEnd'] < len(map):
                            target_seq = query_seq[map[pos['hmmStart']]: map[pos['hmmEnd']] + 1]
                        else:
                            target_seq = ''

                        if re.search('\\A' + condition + '\\Z', target_seq):
                            # we have a pass
                            pass_count += 1

                            # expand possible Nter / Cter positions to seq_from / seq_to
                            if rule['Groups'][grp][pos_num]['start'] == 'Nter':
                                rule['Groups'][grp][pos_num]['start'] = seq_from
                            if rule['Groups'][grp][pos_num]['end'] == 'Cter':
                                rule['Groups'][grp][pos_num]['end'] = seq_to

                    if len(rule['Groups'][grp]) == pass_count:
                        # a group passes only if the whole group is a pass
                        rule_sites.extend(rule['Groups'][grp])

                if rule_sites:
                    domHit = {
                        'score': dom_score,
                        'evalue': dom_evalue,
                        'hmmStart': hmm_from,
                        'hmmEnd': hmm_to,
                        'hmmLength': qlen,
                        'hmmAlign': hmm_align,
                        'start': seq_from,
                        'end': seq_to,
                        'seqAlign': seq_align,
                        'representative': '',
                        'envelopeStart': 1,  # expected result always returns 1
                        'envelopeEnd': 2,  # expected result always returns 2
                        'sites': rule_sites,
                        'scope': rule['Scope'],
                    }

                    if model_id in matches[sequence_id]:
                        matches[sequence_id][model_id]["locations"].append(domHit)
                    else:
                        match_info = {
                            "accession": model_id,
                            "name": model_id,
                            "member_db": member_db,
                            "version": version,
                            "model-ac": model_id.split(":")[0].split(".")[0],
                            "locations": []
                        }
                        matches[sequence_id][model_id] = match_info
                        matches[sequence_id][model_id]["locations"].append(domHit)

    return matches



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

    path_segments = matches.split("/")[-1].split("._.")
    version = path_segments[0]
    member_db = path_segments[1]

    with open(rules_path) as rulesfile:
        rules_hash = json.load(rulesfile)

    matches = matches2rules(matches, rules_hash, member_db, version)

    print(json.dumps(matches, indent=4))
