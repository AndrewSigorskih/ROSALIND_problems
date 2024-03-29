#include "suffix_tree.hpp"

const char* INFILENAME = "rosalind_ba9d.txt";
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


    std::string s1;
    ist >> s1;

    suffix_tree tree(s1);
    tree.lrep(ost);

    return 0;
}