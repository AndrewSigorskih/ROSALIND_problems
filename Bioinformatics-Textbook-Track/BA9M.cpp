#include "bwt.hpp"
//g++ -o ba9m bwt.cpp BA9M.cpp
const char* INFILENAME = "rosalind_ba9m.txt";
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

    std::string transform, pat;
    std::getline(ist, transform);
    BetterBWMatcher matcher(transform);

    while (ist >> pat)
    {
        ost << matcher.match(pat) << " ";
    }
    ost << '\n';
    return 0;
}