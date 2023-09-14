#pragma once
#include <iostream>
#include <fstream>
#include <map>
#include <string>

struct TrieNode;
using childMap = std::map<char, TrieNode*>;

struct TrieNode
{
    TrieNode(uint nodeid) : id(nodeid), isEndOfWord(false) {}
    void printNodes(std::ofstream&);

    uint id;
    bool isEndOfWord;
    childMap children;
};

class Trie
{
public:
    Trie() : nodesNum(0) { root = getNode(); }
    ~Trie();
    void insert(std::string);
    void printTrie(std::ofstream&);
    void TrieMatching(const std::string&, std::ofstream&);
private:
    uint nodesNum;
    TrieNode* root;
    TrieNode* getNode();
    void clearNodes(TrieNode*);
    bool PrefixTrieMatching(const std::string&);
};
