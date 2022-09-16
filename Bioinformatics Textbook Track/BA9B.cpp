#include <iostream>
#include <fstream>
#include <string>
#include <map>

using std::string;
using std::map;
const char* INFILENAME = "rosalind_ba9b.txt";
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
    bool PrefixTrieMatching(const string&);
public:
    Trie() : nodesNum(0) { root = getNode(); }
    ~Trie();
    void insert(string);
    void printTrie(std::ofstream&);
    void TrieMatching(const string&, std::ofstream&);
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

bool Trie::PrefixTrieMatching(const string& text)
{
    TrieNode* curnode = root;
    for(int i = 0; i < text.size(); ++i)
    {
        if (curnode->isEndOfWord)
            return true;
        if (curnode->children.find(text[i]) == curnode->children.end())
            return false;
        curnode = curnode->children[text[i]];
    }
    return curnode->isEndOfWord;
}

void Trie::TrieMatching(const string& text, std::ofstream& ost)
{
    for (int i = 0; i < text.size(); ++i)
    {
        if (PrefixTrieMatching(text.substr(i)))
            ost << i << " ";
    }
    ost << '\n';
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
    string text, tmp;
    getline(ist, text);
    while (getline(ist, tmp))
        trie.insert(tmp);
    trie.TrieMatching(text, ost);
    return 0;
}