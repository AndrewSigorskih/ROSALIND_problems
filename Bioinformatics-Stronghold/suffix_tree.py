# Ukkonen's algorithm  for O(n) suffix tree construction
# source: greck, https://habr.com/ru/post/681940/

# data structure:
# node =  < suffix_link, parent, edges >
# edges = { first_char => edge }
# edge =  < idxs, target_node >
# idxs = (first, last) # corresponds to label = text[first: last]

from collections import defaultdict
from typing import Optional

class Node:
    def __init__(self,
                 edges: Optional[defaultdict]=None,
                 suff_link: Optional['Node']=None,
                 parent: Optional['Node']=None):
        self.edges: defaultdict[str, Edge] = edges if edges is not None else defaultdict(type(None))
        self.suff_link: Optional['Node'] = suff_link
        self.parent: Optional['Node'] = parent
    
    def __str__(self) -> str:
        if self.parent:
            idxs = [e.idxs for e in self.parent.edges.values() if e.target_node == self][0]
            return f"{self.parent}{self.text[idxs[0]: idxs[1]]}_"
        else:
            return "_"

    def is_terminal(self) -> bool:
        return len(self.edges) == 0

class Edge:
    def __init__(self,
                 idxs: Optional[tuple[int, int]]=None,
                 target_node: Optional[Node]=None):
        self.idxs = idxs
        self.target_node = target_node

class SuffixTree:
    def __init__(self, text):
        self.text = text
        self.joker_edge = Edge(idxs=(0, 1))
        self.joker = Node(edges=defaultdict(lambda: self.joker_edge))
        self.root = Node(suff_link=self.joker)
        self.joker_edge.target_node = self.root
        self.infty = len(self.text)
        self._build_tree(self.root, 0, self.infty)

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self.text})'

    def _build_tree(self, node: Node, n: int, infty: int, skip: int = 0):
        while n < infty:
            c = self.text[n]
            edge = node.edges[c]
            if edge is not None:
                first, last = edge.idxs
                i, n0 = first, n
                if skip > 0:
                    can_skip = min(skip, last - first)
                    i += can_skip
                    n += can_skip
                    skip -= can_skip

                while (
                    i < last and n < infty and
                    (self.text[i] == self.text[n] or edge is self.joker_edge)
                ):
                    i += 1
                    n += 1

                if i == last:  # go to the next node
                    node = edge.target_node
                else:  # splitting edge
                    middle_node = Node(parent=node)
                    middle_node.edges[self.text[i]] = edge
                    node.edges[c] = Edge(idxs=(first, i), target_node=middle_node)
                    edge.idxs = (i, edge.idxs[1])
                    edge.target_node.parent = middle_node
                    middle_node.suff_link = self._build_tree(node.suff_link, n0, n, i - first)
                    node = middle_node
            else:  # no way to go; creating new leaf
                new_leaf = Node(parent=node)
                node.edges[c] = Edge(idxs=(n, self.infty), target_node=new_leaf)
                node = node.suff_link
        return node
    
    def get_edges(self, node=None):
        node = node or self.root
        for c, edge in node.edges.items():
            yield self.text[edge.idxs[0]: edge.idxs[1]]
            yield from self.get_edges(edge.target_node)

    def print_tree(self, node=None, indent=0) -> None:
        Node.text = self.text
        node = node or self.root
        space = "    " * indent
        print(f'{space}ID    : {node}')
        print(f'{space}link  : {node.suff_link}')
        print(f'{space}edges : ')

        for c, edge in node.edges.items():
            print(f'{space}  -{c} {edge.idxs}={self.text[edge.idxs[0]: edge.idxs[1]]}:')
            self.print_tree(edge.target_node, indent + 1)

if __name__ == '__main__':
    text = "CCAAGCTGCTAGAGG$CATGCTGGGCTGGCT%"
    tree = SuffixTree(text)
    tree.print_tree()
