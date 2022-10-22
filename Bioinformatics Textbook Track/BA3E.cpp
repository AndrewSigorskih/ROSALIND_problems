#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>

using std::vector;
using std::string;
using std::map;
using grph = map<string, vector<string>>;
const char* INFILENAME = "rosalind_ba3e.txt";
const char* OUTFILENAME = "out.txt";

void dbru(std::ifstream& ist, grph& graph)
{
    string kmer;
    while (getline(ist, kmer))
    {
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

    grph graph;
    dbru(ist, graph);
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