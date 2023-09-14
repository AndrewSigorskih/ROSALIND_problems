#include <iostream>
#include <fstream>
#include <string>
#include <map>

using std::string;
using std::map;
const char* INFILENAME = "rosalind_ba9a.txt";
const char* OUTFILENAME = "out.txt";

class Trie;
class TrieNode;
using childMap = map<char, TrieNode*>;

class TrieNode
{
    uint id;
    bool isEndOfWord;
    childMap children;
    friend class Trie;
public:
    TrieNode(uint nodeid) : id(nodeid), isEndOfWord(false) {}
    void printNodes(std::ofstream&);
    
};

void TrieNode::printNodes(std::ofstream& ost)
{
    for (auto kv: this->children)
    {
        ost << id << "->" << kv.second->id << ":" << kv.first << '\n';
        if (!kv.second->isEndOfWord)
            kv.second->printNodes(ost);
    }
}

class Trie
{
    uint nodesNum;
    TrieNode* root;
    TrieNode* getNode();
    void clearNodes(TrieNode*);
public:
    Trie() : nodesNum(0) { root = getNode(); }
    ~Trie();
    void insert(string);
    void printTrie(std::ofstream&);
};

Trie::~Trie()
{
    clearNodes(root);
}

void Trie::clearNodes(TrieNode* current)
{
    if (!current->isEndOfWord)
    {
        for (auto kv: current->children)
        {
            clearNodes(kv.second);
        }
    }
    delete current;
    --nodesNum;
}

TrieNode* Trie::getNode()
{
    ++nodesNum;
    TrieNode* node = new TrieNode(nodesNum-1);
    return node;
}

void Trie::printTrie(std::ofstream& ost)
{
    root->printNodes(ost);
}

void Trie::insert(string s)
{
    TrieNode* curnode = root;
    for (auto c: s)
    {
        if (curnode->children.find(c) == curnode->children.end())
        {
            curnode->children[c] = getNode();
        }
        curnode = curnode->children[c];
    }
    curnode->isEndOfWord = true;
}

int main()
{
    std::ifstream ist{INFILENAME};
    if (!ist) 
    {
        std::cerr << "Cannot open input file!\n";
        exit(1);
    }
    std::ofstream ost{OUTFILENAME};
    if (!ost) 
    {
        std::cerr << "Cannot open output file!\n";
        exit(1);
    }
    Trie trie;
    string tmp;
    while (getline(ist, tmp))
        trie.insert(tmp);
    trie.printTrie(ost);
    return 0;
}