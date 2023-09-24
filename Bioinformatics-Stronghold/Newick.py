import re
import sys
from collections import deque
from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Optional, TextIO, Tuple, Union

from Branch import posToBits, bitsToPos, get_all_quartets

NWK_REGEX = r"([^:;,()\s]*)(?:\s*:\s*([\d.]+)\s*)?([,);])|(\S)"
Token = Union[str, int]

# parsing inspired by: https://stackoverflow.com/a/51375562

@dataclass(slots=True)
class Node:
    id: int = field(default=0)
    name: str = field(default_factory=str)
    length: float = field(default=1.0, hash=False)
    parentid: int = field(default=-1)
    children: List['Node'] = field(default_factory=list, repr=False)

    def is_leaf(self) -> bool:
        return not self.children

class Tree:
    def __init__(self, newick: Optional[str] = None):
        if not newick:
            self.root = Node(length=0.0)
            return
        if not newick.endswith(';'):
            raise ValueError(f'Newick string should end with ; symbol: {newick}')
        # actual newick parsing
        tokens = re.finditer(NWK_REGEX, newick)
        self.root = self.__parse_node(tokens)[0]
        self.root.length = 0.0

    def __parse_node(self, tokens: Iterator[re.Match],
                     nextid: int = 0, parentid: int = -1) -> Tuple[Node, Token, Token]:
        currId = nextid
        children = []
        name, length, delim, ch = next(tokens).groups(0)
        if ch == "(":
            while ch in "(,":
                child, ch, nextid = self.__parse_node(tokens, nextid+1, currId)
                children.append(child)
            name, length, delim, ch = next(tokens).groups(0)
        return Node(id=currId, name=name,
                    length=float(length) if length else 1.0,
                    parentid=parentid, children=children), delim, nextid

    def __str__ (self) -> str:
        return f'{self.__to_newick(self.root)};'
    
    def __to_newick(self, node: Node, dist: bool=False) -> str:
        children = (
            child.name + (f':{child.length:.3f}'*dist) if child.is_leaf()
            else self.__to_newick(child, dist) + (f':{child.length:.3f}'*dist)
            for child in node.children
        )
        return f'({",".join(children)})'
    
    def to_newick(self, dist: bool=False) -> str:
        return f'{self.__to_newick(self.root, dist)};'

    def print_tree(self, node=None, indent=0, ostream: TextIO=sys.stdout) -> None:
        node = node or self.root
        space = "    " * indent
        print(f'{space}ID : {node.id}, name: {"_" if not node.name else node.name}', file=ostream)
        if node.children:
            print(f'{space}children : ', file=ostream)
            for child in node.children:
                self.print_tree(child, indent+1, ostream)
    
    def travese_preorder(self, node=None) -> Iterator[Node]:
        node = node or self.root
        yield node
        for child in node.children:
            yield from self.travese_preorder(child)

    def traverse_postorder(self, node=None) -> Iterator[Node]:
        node = node or self.root
        for child in node.children:
            yield from self.traverse_postorder(child)
        yield node
    
    def traverse_levelorder(self, node: Optional[Node]=None) -> Iterator[Node]:
        node = node or self.root
        queue = deque()
        queue.append(node)
        while queue:
            node = queue.popleft()
            yield node
            for child in node.children:
                queue.append(child)
    
    def to_branches(self, leavesOrder: Dict[str, int], unrooted: bool = False) -> Dict[int, int]:
        n = len(leavesOrder)
        branches = dict()
        for node in self.traverse_postorder():
            branch = 0
            if node.is_leaf():
                pos = leavesOrder[node.name]
                branch = posToBits(tuple(1 if i == pos else 0 for i in range(n)))
            else:
                for child in node.children:
                    branch |= branches[child.id][0]
            branches[node.id] = (branch, node.is_leaf())
        if unrooted:
            branches.pop(self.root.id)
        return {
            i : branch
            for i, (branch, isLeaf) in branches.items()
            if not isLeaf
        }
    
    @classmethod
    def rf_distance(cls,
                    t1: 'Tree',
                    t2: 'Tree',
                    leavesOrder: Dict[str, int],
                    unrooted: bool = False) -> int:
        substract = 3 if unrooted else 2
        br_1 = set(t1.to_branches(leavesOrder, unrooted).values())
        br_2 = set(t2.to_branches(leavesOrder, unrooted).values())
        return 2*(len(leavesOrder) - substract) - 2*len(br_1 & br_2)    
    
    @classmethod
    def qrt_distance_naive(cls,
                     t1: 'Tree',
                     t2: 'Tree',
                     leavesOrder: Dict[str, int],
                     unrooted: bool = False) -> int:
        from math import comb
        n = len(leavesOrder)
        qrts_1, qrts_2 = set(), set()
        for branch in t1.to_branches(leavesOrder, unrooted).values():
            qrts_1 |= get_all_quartets(bitsToPos(branch, n))
        print(f'{len(qrts_1)=}')
        for branch in t2.to_branches(leavesOrder, unrooted).values():
            qrts_2 |= get_all_quartets(bitsToPos(branch, n))
        print(f'{len(qrts_1)=}')
        return comb(n, 4) * 2 - len(qrts_1 & qrts_2) * 2

class TreeGraph(Tree):
    def __init__(self, newick: Optional[str] = None):
        super().__init__(newick)
        self.nodes: Dict[int, Node] = {}
        self.leaves: Dict[str, Node] = {}
        self.depths: Dict[int, float] = {} # is not initialized by default
        for node in self.travese_preorder():
            self.nodes[node.id] = node
            # assuming all leaf labels are different
            if (node.is_leaf()):
                #if node.name in self.leaves:
                    #raise ValueError(f'Duplicated leaf names are not allowed for {self.__class__.__name__}: {node.name}')
                self.leaves[node.name] = node

    def __name_to_id(self, name: str) -> int:
        if name in self.leaves:
            return self.leaves[name].id
        # the hard way
        for node in self.traverse_levelorder():
            if node.name == name:
                return node.id
        return -1
    
    def __findLCA(self, root: Node, key1: int, key2: int) -> int:
        if root.id == key1 or root.id == key2:
            return root.id
        if root.is_leaf():
            return -1
        child_answers = [self.__findLCA(child, key1, key2) for child in root.children]
        if sum(x > -1 for x in child_answers) >= 2:
            return root.id
        return next((x for x in child_answers if x > -1), -1)
        

    def findLCA(self, key1: Token, key2: Token) -> int:
        if isinstance(key1, str):
            key1 = self.__name_to_id(key1)
        if isinstance(key2, str):
            key2 = self.__name_to_id(key2)
        if any(x == -1 for x in (key1, key2)):
            return -1
        return self.__findLCA(self.root, key1, key2)
    
    def recalc_depths(self)-> None:
        for node in self.travese_preorder():
            parent_depth = 0.0 if node.parentid == -1 else self.depths[node.parentid]
            self.depths[node.id] = parent_depth + node.length

    def get_distance(self, key1: Token, key2: Token) -> float:
        if not self.depths:
            self.recalc_depths()

        if isinstance(key1, str):
            key1 = self.__name_to_id(key1)
        if isinstance(key2, str):
            key2 = self.__name_to_id(key2)

        lca = self.findLCA(key1, key2)
        if any(x == -1 for x in (key1, key2, lca)):
            return -1
        
        return self.depths[key1] + self.depths[key2] - 2*self.depths[lca]
        

if __name__ == '__main__':
    #nwk = "((A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5,G:0.8)F:0.9);"
    nwk = "((((Aa,aa),(Aa,Aa)),((aa,aa),(aa,AA))),Aa);"
    print(nwk)
    tree = Tree(nwk)
    #print(str(tree))
    #print(tree.to_newick(dist=True))
    #tree.print_tree(ostream=sys.stderr)

    print(*(node.id for node in tree.travese_preorder()))
    print(*(node.id for node in tree.traverse_postorder()))
    print(*(node.id for node in tree.traverse_levelorder()))

    #nwk = "((A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5,G:0.8)F:0.9);"
    #g = TreeGraph(nwk)
    #g.print_tree()
    #print(g.leaves)

    nwk = "((A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5,G:0.8)F:0.9);"
    print(nwk)
    t = TreeGraph(nwk)
    t.print_tree()
    print(t.get_distance('A', 'C'))
