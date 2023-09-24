from suffix_tree import SuffixTree

def main():
    with open("rosalind_suff.txt", "r") as f:
        text = f.readline().strip()
    tree = SuffixTree(text)
    with open("out.txt", "w") as o:
       print(*(x for x in tree.get_edges()), sep='\n', file=o)

if __name__ == "__main__":
    main()
