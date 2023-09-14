#include "suffix_array.hpp"
#include <algorithm>

//g++ -o ba9h suffix_array.cpp BA9H.cpp

const char* INFILENAME = "rosalind_ba9h.txt";
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

    std::string text, pat;
    std::getline(ist, text);
    suffix_array array(text);
    std::vector<size_t> result;

    while(std::getline(ist, pat))
    {
        array.all_matches(pat, result);
    }

    std::sort(result.begin(), result.end());
    for (auto i: result)
        ost << i << " ";
    ost << '\n';
    return 0;
}