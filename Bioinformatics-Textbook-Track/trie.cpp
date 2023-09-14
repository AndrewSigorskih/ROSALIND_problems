#include "trie.hpp"

void TrieNode::printNodes(std::ofstream& ost)
{
    for (auto kv: this->children)
    {
        ost << id << "->" << kv.second->id << ":" << kv.first << '\n';
        if (!kv.second->isEndOfWord)
            kv.second->printNodes(ost);
    }
}

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

void Trie::insert(std::string s)
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

bool Trie::PrefixTrieMatching(const std::string& text)
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

void Trie::TrieMatching(const std::string& text, std::ofstream& ost)
{
    for (int i = 0; i < text.size(); ++i)
    {
        if (PrefixTrieMatching(text.substr(i)))
            ost << i << " ";
    }
    ost << '\n';
}
