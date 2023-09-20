#include "bwt.hpp"
//g++ -o ba9n suffix_array.cpp bwt.cpp BA9N.cpp
const char* INFILENAME = "rosalind_ba9n.txt";
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

    std::set<size_t> result;
    std::string text, pat;
    std::getline(ist, text);
    PatternMatcher matcher(text+"$");

    while (ist >> pat)
    {
        matcher.match(pat, result);
    }
    
    for (const auto val: result)
        ost << val << " ";
    ost << '\n';
    return 0;
}