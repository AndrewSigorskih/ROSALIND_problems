#include "suffix_tree.hpp"
// g++ -o ba9e suffix_tree.cpp BA9E.cpp
const char* INFILENAME = "rosalind_ba9e.txt";
const char* OUTFILENAME = "out.txt";

int main()
{
    std::ifstream ist{INFILENAME};
    if (!ist) 
    {
        std::cout << "Cannot open input file!\n";
        exit(1);
    }

    std::ofstream ost{OUTFILENAME};
    if (!ost) 
    {
        std::cout << "Cannot open output file!\n";
        exit(1);
    }


    std::string s1, s2;
    ist >> s1 >> s2;
    size_t first_size = s1.size();
    s1 = s1 + "%" + s2 + "$";

    suffix_tree tree(s1);
    tree.lcs(ost, first_size);
    
    return 0;
}