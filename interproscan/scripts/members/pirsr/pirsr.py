#!/usr/bin/env python

import json
import re
import sys


def matches2rules(matches, rule, member_db, version):
    with open(matches, "r") as matches:
        matches_info = json.load(matches)

    for sequence_id, domains in matches_info.items():
        for model_id, domain in domains.items():
            rule = rules_hash[model_id]
            for location in domain["locations"]:
                hmm_from = location["hmmStart"]
                try:
                    hmm_align = location["hmm_alignment"]
                except KeyError:
                    hmm_align = ""
                seq_align = location["alignment"]

                map = map_hmm_to_seq(hmm_from, hmm_align, seq_align)

                rule_sites = []

                for grp in rule['Groups'].keys():
                    pass_count = 0

                    pos_num = -1
                    for pos in rule['Groups'][grp]:
                        pos_num += 1
                        condition = pos['condition']

                        condition = re.sub('-', '', condition)
                        condition = re.sub('\\(', '{', condition)
                        condition = re.sub('\\)', '}', condition)
                        condition = re.sub('x', '.', condition)

                        query_seq = re.sub('-', '', seq_align)

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
                        # 'score': dom_score,
                        # 'evalue': dom_evalue,
                        # 'hmmStart': hmm_from,
                        # 'hmmEnd': hmm_to,
                        # 'hmmLength': qlen,
                        # 'hmmAlign': hmm_align,
                        # 'start': seq_from,
                        # 'end': seq_to,
                        'seqAlign': seq_align,
                        'representative': '',
                        'envelopeStart': 1,  # expected result always returns 1
                        'envelopeEnd': 2,  # expected result always returns 2
                        'sites': rule_sites,
                        'scope': rule['Scope'],
                    }

                    try:
                        location["domHit"].append(domHit)
                    except KeyError:
                        location["domHit"] = [domHit]

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

    path_segments = matches.split("/")[-1].split("._.")
    version = path_segments[0]
    member_db = path_segments[1]

    with open(rules_path) as rulesfile:
        rules_hash = json.load(rulesfile)

    matches = matches2rules(matches, rules_hash, member_db, version)

    print(json.dumps(matches, indent=4))
