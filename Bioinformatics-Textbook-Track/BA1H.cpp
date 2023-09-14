#include <iostream>
#include <fstream>
#include <string>
#define MIN(a,b) ((a) < (b) ? (a) : (b))
using std::string;
const char* INFILENAME = "rosalind_ba1h.txt";
const char* OUTFILENAME = "out.txt";

int hammingDistance(const string& s1, const string& s2)
{
    if (s1.size() != s2.size())
        return MIN(s1.size(), s2.size());
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
    string pat, text;
    int d;
    ist >> pat >> text >> d;
    for (int i = 0; i <= (text.size() - pat.size()); ++i)
        if (hammingDistance(text.substr(i, pat.size()), pat) <= d)
            ost << i << " ";
    ost << '\n';
    return 0;
}