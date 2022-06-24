import warnings
warnings.filterwarnings("ignore")
#ete3 spams "is with a literal" warnings on import in high python versions
from ete3 import Tree

def main():
    with open("rosalind_sptd.txt", "r") as f:
        _, s1, s2 = [line.strip("\n") for line in f]
    (rf, max_rf, *junk) = Tree(s1).robinson_foulds(Tree(s2), unrooted_trees=True)
    print(rf)

if __name__ == "__main__":
    main()