#include "trie.hpp"
//g++ -o ba9a trie.cpp BA9A.cpp
const char* INFILENAME = "rosalind_ba9a.txt";
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
    std::string tmp;

    while (getline(ist, tmp))
        trie.insert(tmp);

    trie.printTrie(ost);
    return 0;
}