#include "bwt.hpp"
#include <sstream>
//g++ -o ba9o suffix_array.cpp bwt.cpp BA9O.cpp
const char* INFILENAME = "rosalind_ba9o.txt";
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

    size_t mismNum;
    std::string tmp;
    std::vector<size_t> result;
    std::vector<std::string> patterns;


    std::getline(ist, tmp);
    PatternMatcher matcher(tmp+"$");

    std::getline(ist, tmp);
    std::stringstream sstream(tmp);

    while (sstream >> tmp )
    {
        patterns.push_back(tmp);
    }
    ist >> mismNum;

    for (const auto& pat: patterns)
        matcher.match_approx(pat, result, mismNum);
    
    std::sort(result.begin(), result.end());
    for (const auto val: result)
        ost << val << " ";
    ost << '\n';

    return 0;
}