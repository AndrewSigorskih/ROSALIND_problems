#include "trie.hpp"
//g++ -o ba9b trie.cpp BA9B.cpp
const char* INFILENAME = "rosalind_ba9b.txt";
const char* OUTFILENAME = "out.txt";


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
    std::string text, tmp;
    getline(ist, text);

    while (getline(ist, tmp))
        trie.insert(tmp);

    trie.TrieMatching(text, ost);
    return 0;
}