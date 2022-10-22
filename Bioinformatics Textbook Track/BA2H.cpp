#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using std::vector;
using std::string;
const char* INFILENAME = "rosalind_ba2h.txt";
const char* OUTFILENAME = "out.txt";

int hammingDistance(string const& s1, string const& s2)
{
    int res = 0;
    for (int i = 0; i < s1.size(); ++i)
        if (s1[i] != s2[i])
            ++res;
    return res;
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

    string pat, dna;
    getline(ist, pat);
    int k = pat.size(), result = 0, min_dist;
    while (ist >> dna)
    {
        min_dist = INT32_MAX;
        for (int i = 0; i <= (dna.size()-k); ++i)
        {
            int cur_dist = hammingDistance(pat, dna.substr(i, k));
            if (cur_dist < min_dist)
                min_dist = cur_dist;
        }
        result += min_dist;
    }
    ost << result << '\n';
    return 0;
}