#include "bwt.hpp"

//g++ -o ba9i suffix_array.cpp bwt.cpp BA9I.cpp

const char* INFILENAME = "rosalind_ba9i.txt";
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

    std::string text;
    getline(ist, text);
    suffix_array array(text);

    ost << BWTfromArray(text, array.get_array()) << '\n';

    return 0;
}