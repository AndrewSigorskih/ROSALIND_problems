from dataclasses import dataclass, field
from enum import Enum
from typing import List, Mapping

# python because i am tired of parsing...

class Color(Enum):
    GREY = 0
    RED = 1
    BLUE = 2
    PURPLE = 3

    @classmethod
    def string_to_color(self, color: str):
        if (color == 'red'):
            return self.RED
        elif (color == 'blue'):
            return self.BLUE
        else:
            raise ValueError(f'Unknown color value: {color}')
        
    @classmethod
    def color_to_string(self, color):
        if (color == self.GREY):
            return 'grey'
        elif (color == self.RED):
            return 'red'
        elif (color == self.BLUE):
            return 'blue'
        elif (color == self.PURPLE):
            return 'purple'
        else:
            raise ValueError(f'Unknown color value: {color}')

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
        self.colors[id] = Color.string_to_color(color)

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
            color = Color.color_to_string(self.colors[node.id])
            print(f'{node.id}: {color}', file=file)

def main():

    g = Graph()
    with open('rosalind_ba9p.txt', 'r') as f:
        while (line := f.readline().strip()) != '-':
            node, _, children = line.split()
            g.add_node(node, children)
        while (line := f.readline().strip()):
            node, color = line.split(': ')
            g.set_color(node, color)
    g.tree_coloring()
    with open('out.txt', 'w') as f:
        g.print_colors(f)

if __name__ == '__main__':
    main()
