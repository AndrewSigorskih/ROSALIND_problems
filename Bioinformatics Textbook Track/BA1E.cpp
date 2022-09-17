#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>

using std::string;
using std::map;
using std::set;
const char* INFILENAME = "rosalind_ba1e.txt";
const char* OUTFILENAME = "out.txt";

void findClumps(const string& text, int k, int t, set<string>& result)
{
    map<string, int> counts;
    for (int i = 0; (i+k) <= text.size(); ++i)
        counts[text.substr(i, k)] += 1;
    for (auto const& kv: counts)
        if (kv.second >= t)
            result.insert(kv.first);
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

    string text;
    int k, L, t;
    getline(ist, text);
    ist >> k >> L >> t;
    set<string> result;
    for (int i = 0; (i + L) <= text.size(); ++i)
        findClumps(text.substr(i, L), k, t, result);
    for (string x: result)
        ost << x << " ";
    ost << '\n';
    return 0;
}