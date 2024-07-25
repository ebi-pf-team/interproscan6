import re


class re_matcher(object):
    def __init__(self, matchstring):
        self.matchstring = matchstring

    def match(self, regexp):
        self.rematch = re.match(regexp, self.matchstring)
        return bool(self.rematch)

    def group(self, i):
        return self.rematch.group(i)


def stringify(query_id):
    query_id = re.sub(r'[^\w]', '_', query_id)
    return query_id


def parsehmmsearch(hmmer_out):
    matches = {}
    with open(hmmer_out) as fp:
        score_store = {}
        match_store = {}
        forward_seek = False
        domain_hits = {}
        matchpthr = None
        query_id = None

        line = fp.readline()
        while line:
            m = re_matcher(line)
            if line.startswith('#') or not line.strip():
                line = fp.readline()
                continue
            elif m.match(r'\AQuery:\s+(PIRSR[0-9]+\S+)'):
                hmm_id = m.group(1)
                fp.readline()
                fp.readline()
                fp.readline()
                fp.readline()
                line = fp.readline()

                while line.strip():
                    m = re_matcher(line)
                    if m.match(r'\s+------\sinclusion\sthreshold'):
                        break

                    score_array = line.split()
                    if not '---'  in score_array[1]:
                        score_store[stringify(score_array[8])] = float(score_array[1])
                    else:
                        score_store[stringify(score_array[8])] = float(0.0)

                    line = fp.readline()

            elif m.match(r'\A>>\s+(\S+)'):
                query_id = m.group(1)
                query_id = stringify(query_id)

                if not query_id in domain_hits:
                    domain_hits[query_id] = {}
                domain_hits[query_id][hmm_id] = {}

                store_align = 1
                line = fp.readline()

                m = re_matcher(line)
                if m.match(r'\s+\[No individual domains that satisfy reporting thresholds'):
                    break

                fp.readline()
                line = fp.readline()
                while line.strip():
                    domain_info = line.split()
                    if(len(domain_info) != 16):
                        quit()
                    domain_id = domain_info[0]
                    domain_hits[query_id][hmm_id][domain_id] = {}
                    domain_hits[query_id][hmm_id][domain_id]['score'] = domain_info[2]
                    domain_hits[query_id][hmm_id][domain_id]['bias'] = domain_info[3]
                    domain_hits[query_id][hmm_id][domain_id]['cEvalue'] = domain_info[4]
                    domain_hits[query_id][hmm_id][domain_id]['iEvalue'] = domain_info[5]
                    domain_hits[query_id][hmm_id][domain_id]['hmmstart'] = domain_info[6]
                    domain_hits[query_id][hmm_id][domain_id]['hmmend'] = domain_info[7]
                    domain_hits[query_id][hmm_id][domain_id]['alifrom'] = domain_info[9]
                    domain_hits[query_id][hmm_id][domain_id]['alito'] = domain_info[10]
                    domain_hits[query_id][hmm_id][domain_id]['envfrom'] = domain_info[12]
                    domain_hits[query_id][hmm_id][domain_id]['envto'] = domain_info[13]
                    domain_hits[query_id][hmm_id][domain_id]['acc'] = domain_info[15]
                    domain_hits[query_id][hmm_id][domain_id]['hmmalign'] = ''
                    domain_hits[query_id][hmm_id][domain_id]['matchalign'] = ''
                    line = fp.readline()

            elif m.match(r'\s+=='):
                process_align = True
                while process_align:
                    domain_hit_id = m.matchstring.split()[2]
                    hmmalign = ''
                    matchalign = ''
                    store_alignment = False
                    same_domain = True
                    line = fp.readline()
                    while same_domain:
                        hmmalign_array = line.split()
                        hmmalign = hmmalign +  hmmalign_array[2]
                        line = fp.readline()
                        line = fp.readline()
                        matchalign = matchalign + line.split()[2]
                        line = fp.readline()
                        line = fp.readline()
                        line = fp.readline()
                        exit_condition = 0
                        if line.strip():
                            if line.strip().startswith(hmm_id):
                                store_alignment = False
                                dummy = 1
                                exit_condition = 0
                            elif line.strip().startswith('=='):
                                store_alignment = True
                                same_domain = False
                                m = re_matcher(line)
                                exit_condition = 1
                            elif line.strip().startswith('>>'):
                                store_alignment = True
                                m = re_matcher(line)
                                same_domain = False
                                process_align = False
                                exit_condition = 2
                            else:
                                store_alignment = True
                                exit_condition = 3
                        else:
                            store_alignment = True
                            same_domain = False
                            process_align = False
                            exit_condition = 4

                    if store_alignment:
                        domain_hits[query_id][hmm_id][domain_hit_id]['hmmalign'] = hmmalign
                        domain_hits[query_id][hmm_id][domain_hit_id]['matchalign'] = matchalign
                        forward_seek = True

            elif m.match(r'\A\/\/'):
                score_store = {}
            else:
                pass

            if not forward_seek: # we might have peeked forward
                line = fp.readline()
                forward_seek = False
            else:
                forward_seek = False

    fp.close()
    return domain_hits


def parse(output_file):
    matches = parsehmmsearch(output_file)
    raw_matches = []
    for seq_id in matches:
        hits = matches[seq_id]
        for hmm_id  in hits:
            hit = hits[hmm_id]
            for hit_id in hit:
                start = hit[hit_id]['alifrom']
                end = hit[hit_id]['alito']
                matchalign = hit[hit_id]['matchalign']
                hmmalign = hit[hit_id]['hmmalign']
                hmmfrom = hit[hit_id]['hmmstart']
                hmmend = hit[hit_id]['hmmend']
                domscore = hit[hit_id]['score']
                domevalue = hit[hit_id]['iEvalue']
                hit_info = (seq_id, hmm_id, hmmfrom, hmmend, hmmalign, start, end, matchalign, domscore, domevalue)
                raw_matches.append(hit_info)

    return set(raw_matches)
