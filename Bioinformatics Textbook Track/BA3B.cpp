#include <iostream>
#include <fstream>
#include <string>

using std::string;
const char* INFILENAME = "rosalind_ba3b.txt";
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

    string genome, tmp;
    getline(ist, genome);
    while(getline(ist, tmp))
        genome += tmp[tmp.size()-1];
    ost << genome << '\n';
    return 0;
}