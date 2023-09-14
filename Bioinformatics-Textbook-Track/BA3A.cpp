#include <iostream>
#include <fstream>
#include <string>

using std::string;
const char* INFILENAME = "rosalind_ba3a.txt";
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

    int k;
    string text;
    ist >> k >> text;

    for (int i = 0; i < text.size()-k+1; ++i)
        ost << text.substr(i, k) << '\n';
    return 0;
}