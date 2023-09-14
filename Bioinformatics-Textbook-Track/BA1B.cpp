#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>

using std::string;
using std::vector;
using Mymap = std::unordered_map<string, int>;
const char* INFILENAME = "rosalind_ba1b.txt";
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

    string s;
    int k;
    ist >> s >> k;
    Mymap m;
    for (int i = 0; i < s.size()-k; i++)
    {
        m[s.substr(i, k)]++;
    }

    int maxVal = 0;
    for(auto kv : m)
        if(kv.second > maxVal)
            maxVal = kv.second;
    for(auto kv : m)
        if(kv.second == maxVal)
            ost << kv.first << " ";
    ost <<'\n';
    return 0;
}