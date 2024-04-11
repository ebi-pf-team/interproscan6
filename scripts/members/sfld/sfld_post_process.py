import os
import sys
import getopt
import re
from collections import namedtuple

SFLD_POSTPROCESSOR_VERSION = "1.0"

# Named tuples for storing data
HMMERDom = namedtuple('HMMERDom', ['query_ac', 'model_ac', 'seq_evalue', 'seq_score', 'seq_bias', 'dom_cevalue', 'dom_ievalue', 'dom_score', 'dom_bias', 'hmm_start', 'hmm_end', 'ali_start', 'ali_end', 'env_start', 'env_end', 'accuracy'])
SiteMatch = namedtuple('SiteMatch', ['query_ac', 'model_ac', 'n_match_lines', 'match_lines'])
NoHit = namedtuple('NoHit', ['query_ac', 'model_ac'])
Family = namedtuple('Family', ['name', 'n_sites', 'n_features', 'site_pos', 'site_desc', 'site_residue'])
SeqMatches = namedtuple('SeqMatches', ['has_matched', 'residue_matches', 'residue_match_coords'])


def main(argv):
    only_matches = 0
    hmmer_out = None
    output = None
    dom_file = None
    alignments = None
    site_info = None
    format = None

    try:
        opts, args = getopt.getopt(argv, "mp:a:f:d:O:s:o:hv", ["help", "version", "format=", "dom=", "alignments=", "hmmer-out=", "site-info=", "output="])
    except getopt.GetoptError:
        show_help(sys.argv[0])
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-m",):
            only_matches = 1
        elif opt in ("-O",):
            hmmer_out = arg
        elif opt in ("-f",):
            format = arg
        elif opt in ("-o",):
            output = arg
        elif opt in ("-a",):
            alignments = arg
        elif opt in ("-s",):
            site_info = arg
        elif opt in ("-d",):
            dom_file = arg
        elif opt in ("-h", "--help"):
            show_help(sys.argv[0])
            sys.exit()
        elif opt in ("-v", "--version"):
            print(SFLD_POSTPROCESSOR_VERSION)
            sys.exit()

    if not alignments or not hmmer_out or not dom_file or not site_info:
        show_help(sys.argv[0])
        sys.exit(0)

    if output is not None:
        sys.stdout = open(output, "w")

    # Parse HMMER output to find which families have alignments reported
    ali_present, n_families = retrieve_families_with_ali(hmmer_out)

    # Read list of residue-level features
    families = read_site_data(site_info, n_families)

    # Process alignments and check per residue matches
    site_hits, n_site_hits, no_hits, n_no_hits = identify_site_matches(alignments, families, ali_present)

    # Parse domain matches
    dom_hits, n_dom_hits = read_domtblout(dom_file)

    # Remove domain hits to families if there are sequences feature which don't match
    filter_no_hits(dom_hits, n_dom_hits, no_hits, n_no_hits)

    # Having filtered domains, remove the start-end
    strip_dom_se(dom_hits)

    # Output all results
    if format == "text":
        output_dom_sites_as_tab(dom_hits, n_dom_hits, site_hits, n_site_hits)
    else:
        output_dom_sites_by_query(dom_hits, n_dom_hits, site_hits, n_site_hits)

    return 0


def output_dom_sites_by_query(dom_hits, n_dom_hits, site_hits, n_site_hits):
    dom_hits.sort(key=lambda x: (x.query_ac, x.model_ac))
    site_hits.sort(key=lambda x: (x.query_ac, x.model_ac))

    # Skip the matches without accessions
    while n_dom_hits > 0 and not dom_hits[0].query_ac and not dom_hits[0].model_ac:
        n_dom_hits -= 1
        dom_hits = dom_hits[1:]

    h, s = 0, 0
    this_query = ""

    while h < n_dom_hits and s < n_site_hits:
        cmp = cmp_match_pair(dom_hits[h].query_ac, site_hits[s].query_ac, dom_hits[h].model_ac, site_hits[s].model_ac)
        if cmp < 0:  # Next query in alphabetical order is a domain hit
            this_query = dom_hits[h].query_ac
            print(f"Sequence: {this_query}")
            print("Domains:")
            while h < n_dom_hits and dom_hits[h].query_ac == this_query:
                output_dom_hit(dom_hits[h], False)
                h += 1
        elif cmp > 0:  # Next query in alphabetical order is a site hit
            this_query = site_hits[s].query_ac
            print(f"Sequence: {this_query}")
            print("Sites:")
            while s < n_site_hits and site_hits[s].query_ac == this_query:
                output_site_hit(site_hits[s], False)
                s += 1
        else:  # Next query hits to both domains and sites
            this_query = site_hits[s].query_ac
            print(f"Sequence: {this_query}")
            print("Domains:")
            while h < n_dom_hits and dom_hits[h].query_ac == this_query:
                output_dom_hit(dom_hits[h], False)
                h += 1
            print("Sites:")
            while s < n_site_hits and site_hits[s].query_ac == this_query:
                output_site_hit(site_hits[s], False)
                s += 1
        print("// ")

    # Clear up any remaining domain hits...
    while h < n_dom_hits:
        this_query = dom_hits[h].query_ac
        print("// ")
        print(f"Sequence: {this_query}")
        print("Domains:")
        while h < n_dom_hits and dom_hits[h].query_ac == this_query:
            output_dom_hit(dom_hits[h], False)
            h += 1

    # ... or site hits
    while s < n_site_hits:
        this_query = site_hits[s].query_ac
        print("// ")
        print(f"Sequence: {this_query}")
        print("Sites:")
        while s < n_site_hits and site_hits[s].query_ac == this_query:
            output_site_hit(site_hits[s], False)
            s += 1


def output_dom_sites_as_tab(dom_hits, n_dom_hits, site_hits, n_site_hits):
    dom_hits.sort(key=lambda x: (x.query_ac, x.model_ac))
    site_hits.sort(key=lambda x: (x.query_ac, x.model_ac))

    # Skip the matches without accessions
    while n_dom_hits > 0 and not dom_hits[0].query_ac and not dom_hits[0].model_ac:
        n_dom_hits -= 1
        dom_hits = dom_hits[1:]

    for hit in dom_hits:
        output_dom_hit(hit, True)
    for hit in site_hits:
        output_site_hit(hit, True)


def identify_site_matches(aln_fn, families, n_families, ali_present):
    site_hits = []
    n_site_hits = 0
    no_hits = []
    n_no_hits = 0

    with open(aln_fn, 'r') as msaf:
        msa = None
        nfam = 0
        while True:
            msa = read_msa(msaf)
            if msa is None:
                break
            while ali_present[nfam] == 0:
                nfam += 1
            if families[nfam].n_features > 0:
                get_site_matches(families[nfam], msa, site_hits, n_site_hits, no_hits, n_no_hits)
            nfam += 1

    site_hits.sort(key=lambda x: (x.query_ac, x.model_ac))
    no_hits.sort(key=lambda x: (x.query_ac, x.model_ac))
    return site_hits, n_site_hits, no_hits, n_no_hits


def get_site_matches(family, msa, site_hits, n_site_hits, no_hits, n_no_hits):
    rf_array = rf_to_array(msa)
    n_seq = len(msa.sequences)
    pos_map = [get_seq_pos_map(msa.sequences[i]) for i in range(n_seq)]
    seqs_matched = [0] * n_seq
    matches = [SeqMatches(has_matched=0, residue_matches=None, residue_match_coords=None) for _ in range(n_seq)]

    for f in range(family.n_features):
        for s in range(n_seq):
            matches[s].has_matched = 0
            matches[s].residue_matches = [None] * family.n_sites
            matches[s].residue_match_coords = [0] * family.n_sites

        apos = 0
        rpos = 0
        for fi in range(family.n_sites):
            fpos = family.site_pos[fi]
            while apos < len(rf_array) and rf_array[apos] != fpos:
                apos += 1
            if apos < len(rf_array) and rf_array[apos] == fpos:
                for s in range(n_seq):
                    seq_base = msa.sequences[s][apos]
                    if seq_base == family.site_residue[f][fi]:
                        matches[s].residue_matches[rpos] = seq_base
                        matches[s].residue_match_coords[rpos] = pos_map[s][apos]
                        matches[s].has_matched += 1
                rpos += 1
                apos += 1

        for s in range(n_seq):
            if matches[s].has_matched == family.n_sites:
                seqs_matched[s] += 1
                site_hits.append(SiteMatch(
                    query_ac=extract_query_ac(msa.sequences[s]),
                    model_ac=family.name,
                    n_match_lines=build_site_match_strings(family, matches[s].residue_match_coords, matches[s].residue_matches)
                ))
                n_site_hits += 1
            for j in range(family.n_sites):
                free(matches[s].residue_match_coords[j])
                free(matches[s].residue_matches[j])

    for s in range(n_seq):
        if seqs_matched[s] == 0:
            no_hits.append(NoHit(query_ac=extract_query_ac(msa.sequences[s]), model_ac=family.name))
            n_no_hits += 1

    return site_hits, n_site_hits, no_hits, n_no_hits


def build_site_match_strings(family, coords, residues):
    desc = None
    lines = []
    tmp = ""

    while any(c > 0 for c in coords):
        found = False
        for s in range(family.n_sites):
            if coords[s] > 0:
                found = True
                if not desc or desc != family.site_desc[s]:
                    if tmp:
                        lines.append(tmp[:-1] + " " + desc)
                    desc = family.site_desc[s]
                    tmp = ""
                tmp += f"{residues[s]}{coords[s]},"
                coords[s] = 0

        if not found:
            break

    if tmp:
        lines.append(tmp[:-1] + " " + desc)

    return len(lines), lines


def rf_to_array(msa):
    rf_array = [0] * msa.alignment_length
    c = 0
    for i in range(msa.alignment_length):
        if msa.reference_info[i] == 'x':
            rf_array[i] = c
            c += 1
        else:
            rf_array[i] = -1
    return rf_array


def get_seq_pos_map(seq):
    pos_map = [0] * len(seq)
    pos = 0
    for i, c in enumerate(seq):
        if c.isalpha():
            pos_map[i] = pos
            pos += 1
        else:
            pos_map[i] = -1
    return pos_map


def get_options_post(argv):
    long_options = [
        ("help", no_argument, None, 'h'),
        ("version", no_argument, None, 'v'),
        ("format", required_argument, None, 'f'),
        ("dom", required_argument, None, 'd'),
        ("alignments", required_argument, None, 'a'),
        ("hmmer-out", required_argument, None, 'O'),
        ("site-info", required_argument, None, 's'),
        ("output", required_argument, None, 'o'),
    ]

    try:
        opts, args = getopt.getopt(argv, "mp:a:f:d:O:s:o:hv", [o[0] for o in long_options])
    except getopt.GetoptError:
        show_help(sys.argv[0])
        sys.exit(2)

    only_matches = 0
    hmmer_out = None
    output = None
    dom_file = None
    alignments = None
    site_info = None
    format = None

    for opt, arg in opts:
        if opt == "-m":
            only_matches = 1
        elif opt == "-O":
            hmmer_out = arg
        elif opt == "-f":
            format = arg
        elif opt == "-o":
            output = arg
        elif opt == "-a":
            alignments = arg
        elif opt == "-s":
            site_info = arg
        elif opt == "-d":
            dom_file = arg
        elif opt in ("-h", "--help"):
            show_help(sys.argv[0])
            sys.exit()
        elif opt in ("-v", "--version"):
            print(SFLD_POSTPROCESSOR_VERSION)
            sys.exit()

    return only_matches, hmmer_out, output, dom_file, alignments, site_info, format


def get_start_from_nse(s):
    p = s.rindex('/')
    return int(s[p+1:])


def parse_hmmer_dom(line):
    fields = line.split()
    ac = fields[0]
    pd = HMMERDom(
        query_ac=f"{ac}/{int(fields[4].split('/')[1].split('-')[0])}",
        model_ac=fields[3],
        seq_evalue=float(fields[5]),
        seq_score=float(fields[6]),
        seq_bias=float(fields[7]),
        dom_cevalue=float(fields[10]),
        dom_ievalue=float(fields[11]),
        dom_score=float(fields[12]),
        dom_bias=float(fields[13]),
        hmm_start=int(fields[14]),
        hmm_end=int(fields[15]),
        ali_start=int(fields[16]),
        ali_end=int(fields[17]),
        env_start=int(fields[18]),
        env_end=int(fields[19]),
        accuracy=float(fields[20])
    )
    return pd


def read_domtblout(fn):
    dom_hits = []
    n_hits = 0

    with open(fn, 'r') as fp:
        # Skip the header
        while True:
            line = fp.readline()
            if not line.startswith('#'):
                break

        while line:
            dom_hits.append(parse_hmmer_dom(line))
            n_hits += 1
            line = fp.readline()

        if not line.endswith('[ok]'):
            print(f"Failed to see [ok] in file '{fn}': hmmsearch failed?")
            sys.exit(1)

    dom_hits.sort(key=lambda x: (x.query_ac, x.model_ac))
    return dom_hits, n_hits


def filter_no_hits(dom_hits, n_dom_hits, no_hits, n_no_hits):
    n, h = 0, 0
    while h < n_dom_hits:
        cmp = cmp_match_pair(dom_hits[h].query_ac, no_hits[n].query_ac, dom_hits[h].model_ac, no_hits[n].model_ac)
        if cmp > 0:
            n += 1
        elif cmp == 0:
            dom_hits[h] = dom_hits[h]._replace(query_ac='', model_ac='')
            h += 1
            n += 1
        else:
            h += 1


def strip_dom_se(dom_hits):
    for hit in dom_hits:
        if hit.query_ac:
            hit = hit._replace(query_ac=hit.query_ac.split('/')[0])


def retrieve_families_with_ali(fn):
    ali_present = []
    n_families = 0

    with open(fn, 'r') as fp:
        line = fp.readline()
        if not line.startswith('# hmmsearch'):
            print("This does not look like a hmmer output file; expecting '# hmmsearch'")
            sys.exit(1)

        while True:
            line = fp.readline()
            if not line:
                break
            if line.startswith('Accession:'):
                n_families += 1
                ali_present.append(1)
            elif line.startswith('# No hits satisfy inclusion thresholds; no alignment saved'):
                ali_present[-1] = 0

        if not line.endswith('[ok]'):
            print(f"Failed to see [ok] in file '{fn}': hmmsearch failed?")
            sys.exit(1)

    return ali_present, n_families


def read_site_data(fn, n_families):
    families = []

    with open(fn, 'r') as fp:
        line = fp.readline()
        if line != '## MSA feature annotation file\n':
            print("Error reading residue annotation file, incorrect format:")
            print(line)
            sys.exit(1)

        line = fp.readline()
        if not line.startswith('# Format version:'):
            print("Error reading residue annotation file, incorrect format:")
            print(line)
            sys.exit(1)

        if line != '# Format version: 1.1\n':
            print("Error reading residue annotation file: expecting version 1.1")
            sys.exit(1)

        while True:
            line = fp.readline()
            if not line:
                break
            if line.startswith('ACC'):
                fields = line.split()
                name = fields[1]
                n_sites = int(fields[2])
                n_features = int(fields[3])
                family = Family(
                    name=name,
                    n_sites=n_sites,
                    n_features=n_features,
                    site_pos=[0] * n_sites,
                    site_desc=[None] * n_sites,
                    site_residue=[None] * n_features
                )
                for i in range(n_sites):
                    line = fp.readline()
                    if not line.startswith('SITE'):
                        print("Error reading residue annotation file:")
                        print(line)
                        sys.exit(1)
                    fields = line.split()
                    family.site_pos[i] = int(fields[1]) - 1
                    family.site_desc[i] = ' '.join(fields[2:])
                for i in range(n_features):
                    line = fp.readline()
                    if not line.startswith('FEATURE'):
                        print("Error reading residue annotation file:")
                        print(line)
                        sys.exit(1)
                    family.site_residue[i] = list(line[8:8+n_sites])
                families.append(family)

    return families


def output_site_hit(hit, one_line):
    if one_line:
        print(f"{hit.model_ac}", end="")
        for line in hit.match_lines:
            print(f" {line}", end="")
        print(f" {hit.query_ac}")
    else:
        for line in hit.match_lines:
            print(f"{hit.model_ac} {line}")


def output_dom_hit(hit, w_query):
    print(f"{hit.model_ac}\t{hit.seq_evalue:.3e}\t{hit.seq_score:.3e}\t{hit.seq_bias:.3f}\t{hit.hmm_start}\t{hit.hmm_end}\t{hit.dom_score:.3f}\t{hit.ali_start}\t{hit.ali_end}\t{hit.env_start}\t{hit.env_end}\t{hit.dom_cevalue:.3e}\t{hit.dom_ievalue:.3e}\t{hit.accuracy:.3f}\t{hit.dom_bias:.3f}", end="")
    if w_query:
        print(f"\t{hit.query_ac}")
    else:
        print()


def cmp_match_pair(query_1, query_2, model_1, model_2):
    cmp1 = 0 if query_1 == query_2 else 1 if query_1 > query_2 else -1
    cmp2 = 0 if model_1 == model_2 else 1 if model_1 > model_2 else -1
    return cmp1 if cmp1 != 0 else cmp2


def cmp_nohit_query_model(p1, p2):
    return cmp_match_pair(p1.query_ac, p2.query_ac, p1.model_ac, p2.model_ac)


def cmp_dom_query_model(p1, p2):
    return cmp_match_pair(p1.query_ac, p2.query_ac, p1.model_ac, p2.model_ac)


def cmp_site_query_model(p1, p2):
    return cmp_match_pair(p1.query_ac, p2.query_ac, p1.model_ac, p2.model_ac)


def show_help(progname):
    print("Post-process results of HMMER search on SFLD HMMs")
    print(f"Usage {progname}: options:")
    print("\t--version     | -v         show program version")
    print("\t--alignments  | -a         HMMER alignment file")
    print("\t--dom         | -d         HMMER domtblout file")
    print("\t--hmmer-out   | -O         HMMER output file")
    print("\t--site-info   | -s         SFLD reside annotation file")
    print("\t--format      | -f FORMAT  output text format [not implemented]")
    print("\t--output      | -o FILE    output file (otherwise STDOUT)")
    print()


def extract_query_ac(seq_name):
    p = seq_name.rindex('/')
    return seq_name[:p]


def read_msa(fp):
    MSA = namedtuple('MSA', ['sequences', 'reference_info', 'alignment_length'])
    sequences = []
    reference_info = []
    while True:
        line = fp.readline()
        if not line or line.startswith('//'):
            return None if not sequences else MSA(sequences, reference_info, len(sequences[0]))
        if line.startswith('#=GS'):
            continue
        if line.startswith('#=GR'):
            reference_info.append(line.split()[2])
        else:
            sequences.append(line.strip())


if __name__ == "__main__":
    main(sys.argv[1:])
