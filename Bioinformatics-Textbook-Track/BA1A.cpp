#include <iostream>
#include <fstream>
#include <string>

using std::string;
const char* INFILENAME = "rosalind_ba1a.txt";
const char* OUTFILENAME = "out.txt";

int PatternCount(const string& text, const string& pat)
{   // trivial solution
    if (pat.size() > text.size())
        return -1;
    int res = 0;
    for (int i = 0; i <= (text.size() - pat.size()); ++i)
        if (text.substr(i, pat.size()) == pat)
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

    string text, pat;
    ist >> text >> pat;

    ost << PatternCount(text, pat) <<'\n';
    return 0;
}