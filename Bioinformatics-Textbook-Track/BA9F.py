from collections import defaultdict
from enum import Enum
from typing import Iterator, Optional
# needs suffix_tree.py from "stronghold" section
from suffix_tree import SuffixTree, Node
FIRST_STOP = '$'
SECOND_STOP = '%'

class Color(Enum):
    UNKNOWN = 0
    RED = 1
    BLUE = 2
    PURPLE = 3

    def __str__(self):
        return self.name

class ColoredTree(SuffixTree):
    def __init__(self, text):
        super().__init__(text)
        self.colors: defaultdict[int, Color] = defaultdict(lambda: Color.UNKNOWN)
        self.__color_tree()

    def print_tree(self, node=None, indent=0) -> Node:
        Node.text = self.text
        node = node or self.root
        space = "    " * indent
        print(f'{space}ID    : {node}')
        print(f'{space}Color : {self.colors[id(node)]}, leaf: {node.is_terminal()}')
        print(f'{space}edges : ')
        for c, edge in node.edges.items():
            print(f'{space}  -{c} {edge.idxs}={self.text[edge.idxs[0]: edge.idxs[1]]}:')
            self.print_tree(edge.target_node, indent + 1)

    def __color_tree(self, node: Optional[Node] = None) -> Color:
        node = node or self.root
        if not (node.is_terminal()):
            for edge in node.edges.values():
                child_color = self.__color_tree(edge.target_node)
                curr_color = self.colors[id(node)]
                if curr_color == Color.UNKNOWN:
                    self.colors[id(node)] = child_color
                elif (curr_color == Color.RED and child_color == Color.BLUE) or \
                     (curr_color == Color.BLUE and child_color == Color.RED) or \
                     (child_color == Color.PURPLE):
                    self.colors[id(node)] = Color.PURPLE
        else:  # leaf
            idxs = [e.idxs for e in node.parent.edges.values() if e.target_node == node][0]
            substr = self.text[idxs[0]: idxs[1]]
            if (FIRST_STOP in substr):
                self.colors[id(node)] = Color.RED
            else:
                self.colors[id(node)] = Color.BLUE
        return self.colors[id(node)]
    
    def __shorted_uncommon_finder(self, node: Node, seq: list[str]) -> Iterator[str]:
        if self.colors[id(node)] == Color.RED:
            seq[-1] = seq[-1][0]
            result = "".join(seq)
            if (FIRST_STOP not in result) and (SECOND_STOP not in result):
                yield result
        elif self.colors[id(node)] == Color.PURPLE:
            for edge in node.edges.values():
                substr = self.text[edge.idxs[0]: edge.idxs[1]]
                yield from self.__shorted_uncommon_finder(edge.target_node, seq+[substr])

    def shortest_uncommon(self):
        yield from self.__shorted_uncommon_finder(self.root, [])

def main():
    with open('rosalind_ba9f.txt', 'r') as f:
        s1, s2 = f.readline().strip(), f.readline().strip()
    text = f'{s1}{FIRST_STOP}{s2}{SECOND_STOP}'
    tree = ColoredTree(text)
    with open('out.txt', 'w') as o:
        print(
            min(tree.shortest_uncommon(), key=lambda x: len(x)),
            file=o
        )

if __name__ == '__main__':
    main()
