from collections import defaultdict

class SuffixTree():
    def __init__(self, text: str, edges_list: list):
        self.text = text
        self.edges, self.backtrack = dict(), dict()
        self.parents, self.children = set(), set()
        for edge in edges_list:
            parent, child, ind, l = edge.split()
            ind, l = int(ind), int(l)
            self.parents.add(parent)
            self.children.add(child)
            self.backtrack[child] = parent
            self.edges[child] = self.text[(ind - 1) : (ind + l - 1)]

    def _build_seq(self, node: str) -> str:
        res= ""
        while node in self.backtrack:
            res = self.edges[node] + res
            node = self.backtrack[node]
        return res

    def lrep(self, k: int) -> str:
        # estimate lengths of all paths
        paths = defaultdict(int)
        for node in (self.children - self.parents):
            while node in self.backtrack:
                node = self.backtrack[node]
                paths[node] += 1
        
        candidates = [x for x in paths if paths[x] >= k]
        seqs = [self._build_seq(cand) for cand in candidates]
        return max(seqs, key=len)

def main():
    with open("rosalind_lrep.txt", "r") as f:
        text = f.readline().strip()
        k = int(f.readline())
        edges_list = [line.strip() for line in f]
    
    sufftree = SuffixTree(text, edges_list)

    with open("out.txt", "w") as o:
        print(sufftree.lrep(k), file=o)

if __name__ == "__main__":
    main()