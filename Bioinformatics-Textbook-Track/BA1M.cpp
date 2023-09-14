#include <iostream>
#include <fstream>
#include <string>
#include <algorithm> // reverse

using std::string;
const char* INFILENAME = "rosalind_ba1m.txt";
const char* OUTFILENAME = "out.txt";

const char* ALPHABET = "ACGT";

string NumberToPattern(long number, int len)
{
    string result;
    result.reserve(len);
    for(int i = 0; i < len; ++i)
    {
        int tmp = number % 4;
        result += ALPHABET[tmp];
        number = (number-tmp) / 4;
    }
    std::reverse(result.begin(), result.end());
    return result;
}

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

    long N;
    int k;
    ist >> N >> k;

    ost << NumberToPattern(N, k) << '\n';

    return 0;
}