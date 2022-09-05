#include <iostream>
#include <fstream>
#include <string>
#include <map>

using std::map;
using std::string;
const char* INFILENAME = "rosalind_ba1c.txt";
const char* OUTFILENAME = "out.txt";

const map<char, char> Mymap = {
    {'A', 'T'},
    {'T', 'A'},
    {'C', 'G'},
    {'G', 'C'}
};

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
    std::getline(ist, s);
    for (int i = s.size()-1; i >= 0; --i)
    {
        ost << Mymap.at(s[i]);
    }
    ost << '\n';
}