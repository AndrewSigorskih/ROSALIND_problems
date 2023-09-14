import warnings
warnings.filterwarnings("ignore")
#ete3 spams "is with a literal" warnings on import in high python versions
from ete3 import Tree

def tree_to_nontrivial_branches(tree, leaves):
    res = []
    for node in tree.traverse("levelorder"):
        if node.is_leaf():
            continue
        for child in node.get_children():
            if child.is_leaf():
                continue
            child_branch = child.get_leaf_names()
            res.append("".join(["1" if x in child_branch else "0" for x in leaves]))
    return res

def main():
    t = Tree("rosalind_ctbl.txt")
    names = sorted(t.get_leaf_names())
    res = tree_to_nontrivial_branches(t, names)
    with open("out.txt", "w") as o:
        print(*res, sep="\n", file=o)

if __name__ == "__main__":
    main()