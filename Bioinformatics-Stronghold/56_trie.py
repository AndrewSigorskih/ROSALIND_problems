# https://en.wikipedia.org/wiki/Trie
# Classic Trie data structure with fixed size alphabet
class TrieNode:
    rev_dct = {
        0 : "A", 1 : "G",
        2: "C", 3 : "T" 
    }

    def __init__(self, id: int):
        self.children: list = [None]*4
        self.isEndOfWord: bool = False
        self.id: int = id
    
    def printNodes(self, FILE):
        for ind, curnode in enumerate(self.children):
            if curnode:
                print(f"{self.id} {curnode.id} {self.rev_dct[ind]}", file=FILE)
                if not curnode.isEndOfWord:
                    curnode.printNodes(FILE)

class Trie:
    dct = {
        "A" : 0, "G" : 1,
        "C" : 2, "T" : 3, 
    }
    
    def __init__(self):
        self.nodes_num = 0
        self.root = self.getNode()
    
    def getNode(self):
        self.nodes_num += 1
        return TrieNode(self.nodes_num)

    def _charToIndex(self, ch):
        ind = self.dct.get(ch, None)
        if ind is None:
            raise ValueError(f"Wrong alphabet: {ch} is not in (ATGC)")
        return ind
    
    def insert(self, key):
        curnode = self.root
        length = len(key)
        for level in range(length):
            index = self._charToIndex(key[level])
            if not curnode.children[index]:
                curnode.children[index] = self.getNode()
            curnode = curnode.children[index]
        curnode.isEndOfWord = True
    
    def search(self, key):
        curnode = self.root
        length = len(key)
        for level in range(length):
            index = self._charToIndex(key[level])
            if not curnode.children[index]:
                return False
            curnode = curnode.children[index]
        return curnode.isEndOfWord
    
    def printTrie(self, FILE):
        self.root.printNodes(FILE)

def main():
    trie = Trie()
    with open("rosalind_trie.txt", "r") as f:
        for line in f:
            trie.insert(line.strip("\n"))
    
    with open("out.txt", "w") as o:
        trie.printTrie(o)

if __name__ == "__main__":
    main()