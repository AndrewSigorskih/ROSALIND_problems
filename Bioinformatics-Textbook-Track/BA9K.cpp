#include "bwt.hpp"
//g++ -o ba9k bwt.cpp BA9K.cpp

const char* INFILENAME = "rosalind_ba9k.txt";
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

    size_t ind;
    std::string transform;
    std::getline(ist, transform);
    ist >> ind;

    ost << last2first(transform, ind) << '\n';

    return 0;
}