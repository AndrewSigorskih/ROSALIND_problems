#include "suffix_array.hpp"

//g++ -o ba9g suffix_array.cpp BA9G.cpp

const char* INFILENAME = "rosalind_ba9g.txt";
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

    std::string input;
    getline(ist, input);

    suffix_array array(input);
    array.print_array(ost);
    return 0;
}