import json
import sys

MAX_DOM_BY_GROUP = 20
DOM_OVERLAP_THRESHOLD = 0.3
REPR_DOM_DATABASES = ["PFAM", "CDD", "PROSITE_PROFILES", "SMART", "NCBIFAM"]
dc_status_map = {
    "CONTINUOUS": "S",
    "N_TERMINAL_DISC": "N",
    "C_TERMINAL_DISC": "C",
    "NC_TERMINAL_DISC": "NC"
}


def add_representative_domains(matches_path: str) -> list[dict]:
    """
    Add representative domains field
    :param domains: list of domains
    :return: domains with representative domains field added
    """

    with open(matches_path, "r") as f:
        matches = json.load(f)
    for sequence, matches_info in matches.items():
        domains = []
        for acc, info in matches_info.items():
            member_db = info["member_db"]
            for location in info["locations"]:
                location["representative"] = False
                if member_db not in REPR_DOM_DATABASES:
                    continue

                pos_start = location["start"]
                pos_end = location["end"]
                try:
                    fragments = location["location-fragments"]
                    frags_str = ""
                    for fragment in fragments:
                        dc_status = dc_status_map.get(fragment['dc-status'], fragment['dc-status'])
                        frags_str += f"{fragment['start']}-{fragment['end']}-{dc_status},"
                    frags_str = frags_str[:-1]
                except KeyError:
                    frags_str = f"{pos_start}-{pos_end}-S"
                domains.append({
                    "signature": acc,
                    "start": pos_start,
                    "end": pos_end,
                    "frag": frags_str,
                    "fragments": _get_fragments(pos_start, pos_end, frags_str),
                    "rank": REPR_DOM_DATABASES.index(member_db)
                })

            if domains:
                repr_domains = _select_repr_domains(domains)
                for representative in repr_domains:
                    repr_acc = representative["signature"]
                    for location in matches_info[repr_acc]["locations"]:
                        if location["start"] == representative["start"] and location["end"] == representative["end"]:
                            location["representative"] = True
                            break
    return matches


def _select_repr_domains(domains: list[dict]):
    repr_domains = []

    # Sort by boundaries
    domains.sort(key=lambda d: (d["fragments"][0]["start"],
                                d["fragments"][-1]["end"]))

    # Group overlapping domains together
    domain = domains[0]
    domain["residues"] = _calc_coverage(domain)
    stop = domain["fragments"][-1]["end"]
    group = [domain]
    groups = []

    for domain in domains[1:]:
        domain["residues"] = _calc_coverage(domain)
        start = domain["fragments"][0]["start"]

        if start <= stop:
            group.append(domain)
            stop = max(stop, domain["fragments"][-1]["end"])
        else:
            groups.append(group)
            group = [domain]
            stop = domain["fragments"][-1]["end"]

    groups.append(group)

    # Select representative domain in each group
    for group in groups:
        """
        Only consider the "best" N domains of the group, 
        otherwise the number of possible combinations/sets is too high 
        (if M domains, max number of combinations is `2 ^ M`)
        """
        group = sorted(group,
                       key=lambda d: (-len(d["residues"]), d["rank"])
                       )[:MAX_DOM_BY_GROUP]

        nodes = set(range(len(group)))
        graph = {i: nodes - {i} for i in nodes}

        for i, dom_a in enumerate(group):
            for j in range(i + 1, len(group)):
                dom_b = group[j]
                if _eval_overlap(dom_a, dom_b, DOM_OVERLAP_THRESHOLD):
                    graph[i].remove(j)
                    graph[j].remove(i)

        # Find possible domains combinations
        subgroups = _resolve_domains(graph)

        # Find the best combination
        max_coverage = 0
        max_pfams = 0
        best_subgroup = None
        for subgroup in subgroups:
            coverage = set()
            pfams = 0
            _subgroup = []

            for i in subgroup:
                domain = group[i]
                coverage |= domain["residues"]
                if domain["rank"] == 0:
                    pfams += 1

                _subgroup.append(domain)

            coverage = len(coverage)
            if coverage < max_coverage:
                continue
            elif coverage > max_coverage or pfams > max_pfams:
                max_coverage = coverage
                max_pfams = pfams
                best_subgroup = _subgroup

        # Flag selected representative domains
        for domain in best_subgroup:
            repr_domains.append(domain)
    return repr_domains


def _calc_coverage(domain: dict) -> set[int]:
    residues = set()
    for f in domain["fragments"]:
        residues |= set(range(f["start"], f["end"] + 1))

    return residues


def _resolve_domains(graph: dict[int, set[int]]) -> list[set[int]]:
    def is_valid(candidate: list[int]) -> bool:
        for node_a in candidate:
            for node_b in candidate:
                if node_a != node_b and node_a not in graph[node_b]:
                    return False

        return True

    def make_sets(current_set: list[int], remaining_nodes: list[int]):
        if is_valid(current_set):
            if not remaining_nodes:
                all_sets.append(set(current_set))
                return True
        else:
            return False

        current_node = remaining_nodes[0]
        remaining_nodes = remaining_nodes[1:]

        # Explore two possibilities at each step of the recursion
        # 1) current node is added to the set under consideration
        make_sets(current_set + [current_node], remaining_nodes)
        # 2) current node is not added to the set
        make_sets(current_set, remaining_nodes)

    all_sets = []
    make_sets([], list(graph.keys()))
    return all_sets


def _eval_overlap(dom_a: dict, dom_b: dict, threshold: float) -> bool:
    overlap = len(dom_a["residues"] & dom_b["residues"])
    return overlap and overlap / min(len(dom_a["residues"]),
                                     len(dom_b["residues"])) >= threshold


def _get_fragments(pos_start: int, pos_end: int, fragments: str) -> list[dict]:
    if fragments:
        result = []
        for frag in fragments.split(','):
            # Format: START-END-STATUS
            s, e, t = frag.split('-')
            result.append({
                "start": int(s),
                "end": int(e),
                "dc-status": t
            })

        result.sort(key=lambda x: (x["start"], x["end"]))
    else:
        result = [{
            "start": pos_start,
            "end": pos_end,
            "dc-status": "S"  # Continuous
        }]

    return result


def main():
    args = sys.argv[1:]
    matches_path = args[0]
    matches_with_repr = add_representative_domains(matches_path)
    print(json.dumps(matches_with_repr, indent=2))


if __name__ == '__main__':
    main()
