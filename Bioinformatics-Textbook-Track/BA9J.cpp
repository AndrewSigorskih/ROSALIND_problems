#include "bwt.hpp"
//g++ -o ba9j bwt.cpp BA9J.cpp

const char* INFILENAME = "rosalind_ba9j.txt";
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

    std::string transform;
    std::getline(ist, transform);
    ost << inverse_bwt(transform) << '\n';

    return 0;
}