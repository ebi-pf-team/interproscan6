import logging
import os
import sys
import getopt

from collections import namedtuple
from pathlib import Path

from sfld_post_process_utilities import (
    Family,
    SiteMatch,
    SFLD_POSTPROCESSOR_VERSION,
    build_parser,
)


logger = logging.getLogger()


def main():
    ali_present = None

    parser = build_parser()
    args = parser.parse_args()

    if args.version:
        print(SFLD_POSTPROCESSOR_VERSION)
        sys.exit(0)

    if None in (args.alignment, args.dom, args.hmmer_out, args.site_info):
        logger.error(
            (
                "The alignment file (--alignment), domtblout (--dom) file, "
                " HMMER out file (--hmmer_out), and site annotation file (--site_info)"
                " must all be defined\n"
                "Terminating program"
            )
        )
        sys.exit(1)

    if not check_hmmsearch_status(args.hmmer_out):
        print(
            f"Failed to see [ok] in file '{args.hmmer_out}': hmmsearch failed?"
        )
        sys.exit(1)

    # Parse HMMER output to find which interpro families have alignments reported
    # {acc: 1 = alignment, 0 = no alignemnt}
    ali_present = retrieve_families_with_ali(args.hmmer_out)

    # Read list of residue-level features
    interpro_families = read_site_data(args.site_info)

    # Process alignments and check per residue matches
    site_hits, n_site_hits, no_hits, n_no_hits = identify_site_matches(
        interpro_families,
        ali_present
    )

    # # print(site_hits)
    # # print(n_site_hits)
    # # print(no_hits)
    # # print(n_no_hits)

    # # # Parse domain matches
    # # dom_hits, n_dom_hits = read_domtblout(dom_file)

    # # # Remove domain hits to families if there are sequences feature which don't match
    # # filter_no_hits(dom_hits, n_dom_hits, no_hits, n_no_hits)

    # # # Having filtered domains, remove the start-end
    # # strip_dom_se(dom_hits)

    # # # Output all results
    # # if format == "text":
    # #     output_dom_sites_as_tab(dom_hits, n_dom_hits, site_hits, n_site_hits)
    # # else:
    # #     output_dom_sites_by_query(dom_hits, n_dom_hits, site_hits, n_site_hits)

    # # return 0


def check_hmmsearch_status(hmmer_out: Path) -> bool:
    """Check hmmsearch was completed ok"""
    with open(hmmer_out, 'r') as f:
        lines = f.readlines()
        if len(lines) < 1:
            return False
        last_line = lines[-1].strip()
        if last_line == "[ok]":
            return True
        if len(lines) >= 2:
            second_last_line = lines[-2].strip()
            if second_last_line == "[ok]":
                return True
        return False


def retrieve_families_with_ali(hmmer_out: Path) -> dict[str, dict[str, dict[str, list[str]]]]:
    """Parse HMMER output to find which interpro families have alignments reported

    :param hmmer_out: path to HMMER output (.out) file

    Marks if alignment is present for accession -- retrieves lies containing alignemnt
    sequence for the interpro family and the query protein sequence into a list
    E.g. {'SFLDS00029':
        {'family': {'query_acc1': ['2 vvtrgC 7']}, 'query': {'query_acc1': ['144 VVcarGC 150']}}}
    """
    ali_present = {}
    acc, query_acc = 'accession_placeholder', 'query_acc_placeholder'

    with open(hmmer_out, 'r') as fp:
        line = fp.readline()
        if not line.startswith('# hmmsearch'):
            print(
                "This does not look like a hmmer output file; expecting '# hmmsearch' at start")
            sys.exit(1)

        while True:
            line = fp.readline()
            if not line:
                break
            if line.startswith('Accession:'):
                acc = line.split(" ")[-1].strip("\n")
                # + one key per query protein with a sig hit#
                ali_present[acc] = {'family': {}, 'query': {}}

            elif line.startswith('>> '):
                query_acc = line.split()[1]

            elif line.strip().startswith(acc):
                try:
                    ali_present[acc]['family'][query_acc].append(
                        " ".join(line.strip().split()[1:]))
                except KeyError:
                    ali_present[acc]['family'][query_acc] = [
                        " ".join(line.strip().split()[1:])]

            elif line.strip().startswith(query_acc):
                try:
                    ali_present[acc]['query'][query_acc].append(
                        " ".join(line.strip().split()[1:]))
                except KeyError:
                    ali_present[acc]['query'][query_acc] = [
                        " ".join(line.strip().split()[1:])]

    return ali_present


def read_site_data(site_info: Path) -> dict[str, Family]:
    """Read list of residue-level features.

    Returns a list of Family objects, which represent as family
    as described in the site_info file.

    :param site_info: path to SFLD site annotation file

    Example Family:
        NAME: SFLDG01108
        n_sites: 3
        n_features: 1
        site_pos: [37, 41, 44]
        site_desc: ['Binds [4Fe-4S]-AdoMet cluster', 'Binds [4Fe-4S]-AdoMet cluster', 'Binds [4Fe-4S]-AdoMet cluster']
        site_residue: ['C', 'C', 'C']
        sequence: [<start pos> <seq> <end pos>]
    """
    families = {}

    with open(site_info, 'r') as fp:
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
                fields = line.split()  # e.g. ['ACC', 'SFLDS00001', '3', '8']
                name = fields[1]
                n_sites = int(fields[2])
                n_features = int(fields[3])

                families[name] = Family()
                families[name].name = fields[1]
                families[name].n_sites = n_sites
                families[name].n_features = n_features
                families[name].site_pos = [0] * families[name].n_sites
                families[name].site_desc = [None] * families[name].n_sites
                families[name].site_residue = [None] * families[name].n_sites

                for i in range(n_sites):
                    line = fp.readline()
                    if not line.startswith('SITE'):
                        logger.error(
                            "Error reading residue annotation file: %s", line)
                        sys.exit(1)
                    fields = line.split()
                    families[name].site_pos[i] = int(fields[1]) - 1
                    families[name].site_desc[i] = ' '.join(fields[2:])
                for i in range(n_features):
                    line = fp.readline()
                    if not line.startswith('FEATURE'):
                        logger.error(
                            "Error reading residue annotation file: %s", line)
                        sys.exit(1)
                    families[name].site_residue = list(line[8:8+n_sites])

    return families


def identify_site_matches(
    interpro_families: dict[str, Family],
    ali_present: dict[str, list[str]]
) -> list[SiteMatch]:
    """
    Screen the domain hits retrieved from the hmmer.out file, and see if any of the sites
    retrieved from the site annotation file are conserved in the domain hits --> SiteMatch

    :param interpro_families: dict of families from the site annotation file
    :param ali_present: marks for each protein if alignment is present in the hmmer.out file
    """
    matches = []

    for family in interpro_families:
        # no alignments = no significant hits found
        if len(ali_present[family]['family']) == 0:
            continue

        # check if domain hit contains site hits
        for i, (site_pos, site_residue) in enumerate(
            zip(interpro_families[family].site_pos, interpro_families[family].site_residue)
        ):
            for query_acc in ali_present[family]['family']:  # query protein ID
                for j, algn_line in enumerate(ali_present[family]['family'][query_acc]):
                    algn_start, algn_end = int(algn_line.split()[0]), int(algn_line.split()[-1])

                    if algn_start <= site_pos <= algn_end:
                        # check if interpro site in alignment line
                        site_index = site_pos - algn_start + 1
                        query_algn = ali_present[family]['query'][query_acc][j]
                        query_residue = query_algn.split()[1][site_index]

                        if query_residue.lower() == site_residue.lower():
                            hit = SiteMatch()
                            hit.query_ac = query_acc
                            hit.model_ac = family
                            hit.site_pos = site_index + int(query_algn.split()[0])
                            hit.site_desc = interpro_families[family].site_desc[i]
                            hit.site_residue = query_residue
                            matches.append(hit)

    return matches


def get_site_matches(
    family: Family,
    msa: namedtuple,
    site_hits: list,
    n_site_hits: int,
    no_hits: list,
    n_no_hits: int
) -> tuple[list, int, list, int]:
    """
    Get per residue matches using alignments for this family
    """
    rf_array = rf_to_array(msa)
    n_seq = len(msa.sequences)
    pos_map = [get_seq_pos_map(msa.sequences[i]) for i in range(n_seq)]
    seqs_matched = [0] * n_seq
    matches = [SeqMatches() for _ in range(n_seq)]

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
                site_match = SiteMatch()
                site_match.query_ac = extract_query_ac(msa.sequences[s])
                site_match.model_ac = family.name
                site_match.n_match_lines = build_site_match_strings(
                    family,
                    matches[s].residue_match_coords,
                    matches[s].residue_matches
                )
                site_hits.append(site_match)
                n_site_hits += 1

    for s in range(n_seq):
        if seqs_matched[s] == 0:
            no_hit = NoHit()
            no_hit.query_ac = extract_query_ac(msa.sequences[s])
            no_hit.model_ac = family.name
            no_hits.append(no_hit)
            n_no_hits += 1

    return site_hits, n_site_hits, no_hits, n_no_hits


def rf_to_array(msa: namedtuple) -> list[int]:
    rf_array = [0] * msa.alignment_length
    c = 0
    for i in range(msa.alignment_length):
        if msa.reference_info[i] == 'x':
            rf_array[i] = c
            c += 1
        else:
            rf_array[i] = -1
    return rf_array


def extract_query_ac(seq_name: str) -> str:
    p = seq_name.rindex('/')
    return seq_name[:p]


def build_site_match_strings(
    family: Family,
    coords,
    residues
) -> tuple[int, list[str]]:
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
        cmp = cmp_match_pair(
            dom_hits[h].query_ac, site_hits[s].query_ac, dom_hits[h].model_ac, site_hits[s].model_ac)
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
        opts, args = getopt.getopt(argv, "mp:a:f:d:O:s:o:hv", [
                                   o[0] for o in long_options])
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


def get_start_from_nse(s: str) -> int:
    """Retrieves the start position of ???"""
    p = s.rindex('/')
    return int(s[p+1:])


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


def filter_no_hits(dom_hits, n_dom_hits, no_hits, n_no_hits):
    n, h = 0, 0
    while h < n_dom_hits:
        cmp = cmp_match_pair(
            dom_hits[h].query_ac, no_hits[n].query_ac, dom_hits[h].model_ac, no_hits[n].model_ac)
        if cmp > 0:
            n += 1
        elif cmp == 0:
            dom_hits[h] = dom_hits[h]._replace(query_ac='', model_ac='')
            h += 1
            n += 1
        else:
            h += 1


def cmp_match_pair(query_1, query_2, model_1, model_2):
    cmp1 = 0 if query_1 == query_2 else 1 if query_1 > query_2 else -1
    cmp2 = 0 if model_1 == model_2 else 1 if model_1 > model_2 else -1
    return cmp1 if cmp1 != 0 else cmp2


def strip_dom_se(dom_hits):
    for hit in dom_hits:
        if hit.query_ac:
            hit = hit._replace(query_ac=hit.query_ac.split('/')[0])


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


if __name__ == "__main__":
    main()
