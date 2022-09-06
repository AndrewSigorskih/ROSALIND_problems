#include <iostream>
#include <fstream>
#include <string>
#include <map>

using std::map;
using std::string;
const char* INFILENAME = "rosalind_ba1l.txt";
const char* OUTFILENAME = "out.txt";

const map<char, int> Mymap = {
    {'A', 0},
    {'C', 1},
    {'G', 2},
    {'T', 3}
};

long PatternToNumber(const string& s)
{
    if (s.size() == 0)
        return -1L;
    long result = 0;
    for (char c: s)
        result = 4 * result + Mymap.at(c);
    return result;
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

    string s;
    getline(ist, s);
    ost << PatternToNumber(s) << '\n';

    return 0;
}