from collections import defaultdict
from typing import List, Optional

# ----------------------------- #
# Construct Suffix Tree from Suffix Array and LCP Array
# https://stackoverflow.com/questions/57502012/how-to-construct-suffix-tree-from-lcp-array-and-suffix-array
# https://en.wikipedia.org/wiki/LCP_array#Suffix_tree_construction
# ----------------------------- #

class Node:
    def __init__(self,
                 depth: Optional[int] = None):
        self._depth = depth
        self.parent = None
        self.parent_edge = None
        self.edges: defaultdict[str, Edge] = defaultdict(type(None))

    @property
    def depth(self) -> int:
        if self._depth is None:
            self.__calc_depth()
        return self._depth
    
    def unset_depth(self):
        self._depth = None
    
    def __calc_depth(self) -> None:
        if self.parent is None:
            self._depth = 0
        else:
            self._depth = self.parent.depth + len(self.parent_edge)


class Edge:
    def __init__(self,
                 idxs: Optional[tuple[int, int]]=None,
                 target_node: Optional[Node]=None):
        self.idxs = idxs
        self.target_node = target_node

    def __len__(self) -> int:
        return self.idxs[1] - self.idxs[0]
    
    def get_label(self, text: str):
        return text[self.idxs[0]: self.idxs[1]]

class Tree:
    def link_nodes(self, 
                   parent: Node, child:Node,
                   start: int, end: int):
        edge = Edge((start, end), child)
        child.parent = parent
        parent.edges[self.text[start]] = edge
        child.parent_edge = edge

    def __init__(self, text: str, sarr: List[int], lcp: List[int]):
        self.root = Node()
        self.text = text
        # make p a child of root with edge label S[0] (the least suffix)
        p = Node()
        edge = Edge((sarr[0], len(text)), p)
        self.link_nodes(self.root, p, sarr[0], len(text))

        prevP = None
        for i in range(1, len(text)):
            l = lcp[i]
            while p.depth > l:
                prevP = p
                p = p.parent
            if p.depth == l:
                # make q a child of p, with edge label S[i][l...]
                q = Node()
                self.link_nodes(p, q, sarr[i]+l, len(text))
                # p := q
                p = q
            else:
                # we need to "split" the edge from p to prevP
                prevDepth = prevP.depth
                # unlink prevP from p
                p.edges.pop(self.text[prevP.parent_edge.idxs[0]])
                prevP.unset_depth()
                # make q a child of p, with edge label S[i-1][(depth of p)...(l - 1)]
                q = Node()
                pDepth = p.depth
                self.link_nodes(p, q, sarr[i-1]+pDepth, sarr[i-1]+l)
                # make prevP a child of q, with edge label S[i-1][l...(prevDepth - 1)]
                self.link_nodes(q, prevP, sarr[i-1]+l, sarr[i-1]+prevDepth)
                # make r a child of q, with edge label S[i][l...]
                r = Node()
                self.link_nodes(q, r, sarr[i]+l, len(text))
                # p := r
                p = r

    def get_edges(self, node=None):
        node = node or self.root
        for _, edge in node.edges.items():
            yield edge.get_label(self.text)
            yield from self.get_edges(edge.target_node)


def main():
    with open('rosalind_ba9r.txt', 'r') as f:
        text = f.readline().strip()
        sarr = [int(x) for x in f.readline().strip().split(', ')]
        lcp = [int(x) for x in f.readline().strip().split(', ')]
    with open('out.txt', 'w') as o:
        print(*Tree(text, sarr, lcp).get_edges(), sep='\n', file=o)

if __name__ == '__main__':
    main()
