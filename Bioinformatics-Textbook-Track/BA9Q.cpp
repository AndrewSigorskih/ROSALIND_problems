#include "suffix_array.hpp"
//g++ -o ba9q suffix_array.cpp BA9Q.cpp
const char* INFILENAME = "rosalind_ba9q.txt";
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

    int K;
    std::string text;
    std::getline(ist, text);
    ist >> K;

    suffix_array sarr(text);
    const std::vector<size_t>& arr = sarr.get_array();

    for (size_t i = 0; i < arr.size(); ++i)
    {
        if (arr[i] % K == 0)
            ost << i << "," << arr[i] << '\n';
    }

    return 0;
}