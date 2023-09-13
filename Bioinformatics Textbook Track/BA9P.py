from dataclasses import dataclass, field
from enum import Enum
from typing import List, Mapping

# python because i am tired of parsing...

class Color(Enum):
    GREY = 0
    RED = 1
    BLUE = 2
    PURPLE = 3

    def __str__(self) -> str:
        return self.name.lower()
    
COLOR_MAPPING = {
    'red': Color.RED,
    'blue': Color.BLUE,
    'purple': Color.PURPLE,
    'grey': Color.GREY
}

def string_to_color(s: str) -> Color:
    if s not in COLOR_MAPPING:
        raise ValueError(f'Unknown color value: {s}')
    return COLOR_MAPPING[s]


@dataclass
class Node:
    id: int
    children: List[int] = field(default_factory=list)

    def isLeaf(self) -> bool:
        return len(self.children) == 0

@dataclass
class Graph:
    nodes: List[Node] = field(default_factory=list)
    colors: Mapping[int, Color] = field(default_factory=dict)   

    def add_node(self, id: str, children: str) -> None:
        id = int(id)
        children = None if children == "{}" else list(map(int, children.split(',')))
        self.nodes.append(
            Node(id) if not children else Node(id, children)
        )
        self.colors[id] = Color.GREY

    def set_color(self, id: str, color: str) -> None:
        id = int(id)
        self.colors[id] = string_to_color(color)

    def tree_coloring(self) -> None:
        grey_nodes = [node.id for node in self.nodes 
                      if self.colors[node.id] == Color.GREY]
        for id in grey_nodes:
            if self.colors[id] != Color.GREY:  # may be set during recursive call
                continue
            self.colors[id] = self.__color_node(id)

    def __color_node(self, id: int) -> Color:
        if self.nodes[id].isLeaf():
            return self.colors[id]
        child_colors = set()
        for child in self.nodes[id].children:
            if self.colors[child] == Color.GREY:
                self.colors[child] = self.__color_node(child)
            child_colors.add(self.colors[child])
        if len(child_colors) >= 2:
            return Color.PURPLE
        else:
            return child_colors.pop()
        
    def print_colors(self, file):
        for node in self.nodes:
            print(f'{node.id}: {self.colors[node.id]}', file=file)

def main():

    graph = Graph()
    with open('rosalind_ba9p.txt', 'r') as f:
        while (line := f.readline().strip()) != '-':
            node, _, children = line.split()
            graph.add_node(node, children)
        while (line := f.readline().strip()):
            node, color = line.split(': ')
            graph.set_color(node, color)
    graph.tree_coloring()
    with open('out.txt', 'w') as f:
        graph.print_colors(f)

if __name__ == '__main__':
    main()
