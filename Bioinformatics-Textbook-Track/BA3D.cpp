#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>

using std::vector;
using std::string;
using std::map;
using grph = map<string, vector<string>>;
const char* INFILENAME = "rosalind_ba3d.txt";
const char* OUTFILENAME = "out.txt";

void dbru(string const& text, int k, grph& graph)
{
    for (int i = 0; i <= (text.size()-k); ++i)
    {
        string kmer = text.substr(i, k);
        graph[kmer.substr(0, kmer.size()-1)].push_back(kmer.substr(1));
    }
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

    int k;
    string text;
    ist >> k >> text;

    grph graph;
    dbru(text, k, graph);
    for (auto kv: graph)
    {
        ost << kv.first << " -> ";
        ost << kv.second[0];
        for (int i = 1; i < kv.second.size(); ++i)
            ost << "," << kv.second[i];
        ost << '\n';
    }
    return 0;
}