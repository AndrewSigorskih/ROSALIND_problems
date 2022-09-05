#include <iostream>
#include <fstream>
#include <string>

using std::string;
const char* INFILENAME = "rosalind_ba1g.txt";
const char* OUTFILENAME = "out.txt";

int hammingDistance(const string& s1, const string& s2)
{
    if (s1.size() != s2.size())
        return -1;
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
    string s1, s2;
    std::getline(ist, s1);
    std::getline(ist, s2);
    ost << hammingDistance(s1, s2) << '\n';
    return 0;
}