import json
import re
import sys

from Bio import Phylo
from Bio.Phylo import NewickIO


def main():
    jplace_file = sys.argv[1]
    newick_file = sys.argv[2]

    with open(jplace_file, "rt") as fh:
        results = json.load(fh)

    tree_string = results["tree"]
    matches = re.findall(r"AN(\d+):\d+\.\d+\{(\d+)\}", tree_string)

    an_label = {}
    for [an, r] in matches:
        an_label["AN" + an] = "R" + r
        an_label["R" + r] = "AN" + an

    newick_string = re.sub(r"(AN\d+)?\:\d+\.\d+{(\d+)}", r"R\g<2>",
                           tree_string)

    newick_string = re.sub(r"AN\d+", r"", newick_string)
    newick_string = re.sub(r"BI\d+", r"", newick_string)
    mytree = Phylo.read(NewickIO.StringIO(newick_string), "newick")
    
    for placement in results["placements"]:
        query_id = placement["n"][0]
        child_ids = []
        ter = []

        for maploc in placement["p"]:
            rloc = "R" + str(maploc[0])
            clade_obj = mytree.find_clades(rloc)

            node = next(clade_obj)
            ter.extend(node.get_terminals())
            comonancestor = mytree.common_ancestor(ter)

            for leaf in comonancestor.get_terminals():
                child_ids.append(an_label[leaf.name])

        newtree = Phylo.read(newick_file, "newick")
        common_an = newtree.common_ancestor(child_ids)
        print(query_id, str(common_an) if common_an else "root")


if __name__ == "__main__":
    main()
